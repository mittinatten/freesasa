#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdlib.h>

#include "classifier.h"
#include "coord.h"
#include "freesasa_internal.h"
#include "pdb.h"

/**
   This file contains the functions that turn lines in PDB files to
   atom records. It's all pretty messy because one needs to keep track
   of when a new chain/new residue is encountered, and skip atoms that
   are duplicates (only first of alt coordinates are used). It is both
   possible to add atoms one by one, or by reading a whole file at
   once, and the user can supply some options about what to do when
   encountering atoms that are not recognized, whether to include
   hydrogrens, hetatm, etc. The current implementation results in a
   rather convoluted logic, that is difficult to maintain, debug and
   extend.

   TODO: Refactor.
*/
#define ATOMS_CHUNK 512
#define RESIDUES_CHUNK 64
#define CHAINS_CHUNK 64
#define CHAIN_LABEL_LENGTH 3

typedef char chain_label_t[CHAIN_LABEL_LENGTH + 1];

struct atom {
    char res_name[PDB_ATOM_RES_NAME_STRL + 1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL + 1];
    char atom_name[PDB_ATOM_NAME_STRL + 1];
    char symbol[PDB_ATOM_SYMBOL_STRL + 1];
    char *line;
    int res_index;
    chain_label_t chain_label;
    freesasa_atom_class the_class;
};

static const struct atom empty_atom = {
    "",                   /* res_name */
    "",                   /* res_number */
    "",                   /* atom_name */
    "",                   /* symbol */
    NULL,                 /* line */
    -1,                   /* res_index */
    "",                   /* chain_label */
    FREESASA_ATOM_UNKNOWN /* the_class */
};

struct atoms {
    int n;
    int n_alloc;
    struct atom **atom;
    double *radius;
};

struct residues {
    int n;
    int n_alloc;
    int *first_atom;
    freesasa_nodearea **reference_area;
};

struct chains {
    int n;
    int n_alloc;
    chain_label_t *labels;
    char *short_labels; /* NULL-terminated string */
    int *first_atom;    /* first atom of each chain */
};

struct freesasa_structure {
    struct atoms atoms;
    struct residues residues;
    struct chains chains;
    char *classifier_name;
    coord_t *xyz;
    int model; /* model number */
    size_t cif_ref;
    void (*release_cif_ref)(size_t);
};

static int
guess_symbol(char *symbol,
             const char *name);

static void
atom_free(struct atom *a)
{
    if (a != NULL) {
        free(a->line);
    }
    free(a);
}

struct atoms
atoms_init(void)
{
    struct atoms atoms;
    atoms.n = 0;
    atoms.n_alloc = 0;
    atoms.atom = NULL;
    atoms.radius = NULL;
    return atoms;
}

/* Allocates memory in chunks, ticks up atoms->n if allocation successful */
static int
atoms_alloc(struct atoms *atoms)
{
    int i, new_size;
    void *aa, *ar;

    assert(atoms);
    assert(atoms->n <= atoms->n_alloc);

    if (atoms->n == atoms->n_alloc) {
        new_size = atoms->n_alloc + ATOMS_CHUNK;
        aa = atoms->atom;
        ar = atoms->radius;

        atoms->atom = realloc(atoms->atom, sizeof(struct atom *) * new_size);
        if (atoms->atom == NULL) {
            atoms->atom = aa;
            return mem_fail();
        }

        for (i = atoms->n_alloc; i < new_size; ++i) {
            atoms->atom[i] = NULL;
        }

        atoms->radius = realloc(atoms->radius, sizeof(double) * new_size);
        if (atoms->radius == NULL) {
            atoms->radius = ar;
            return mem_fail();
        }

        atoms->n_alloc = new_size;
    }
    ++atoms->n;
    return FREESASA_SUCCESS;
}

static void
atoms_dealloc(struct atoms *atoms)
{
    int i;
    struct atom **atom;

    if (atoms) {
        atom = atoms->atom;
        if (atom) {
            for (i = 0; i < atoms->n; ++i)
                if (atom[i]) atom_free(atom[i]);
            free(atom);
        }
        free(atoms->radius);
        *atoms = atoms_init();
    }
}

static struct atom *
atom_new(const char *residue_name,
         const char *residue_number,
         const char *atom_name,
         const char *symbol,
         const chain_label_t chain_label)
{
    struct atom *a = malloc(sizeof(struct atom));

    if (a == NULL) {
        mem_fail();
        free(a);
    } else {
        *a = empty_atom;

        a->line = NULL;
        a->res_index = -1;

        snprintf(a->atom_name, sizeof(a->atom_name), "%s", atom_name);
        snprintf(a->res_name, sizeof(a->res_name), "%s", residue_name);
        snprintf(a->res_number, sizeof(a->res_number), "%s", residue_number);
        snprintf(a->symbol, sizeof(a->symbol), "%s", symbol);
        snprintf(a->chain_label, sizeof(chain_label_t), "%s", chain_label);

        a->the_class = FREESASA_ATOM_UNKNOWN;
    }

    return a;
}

static struct atom *
atom_new_from_line(const char *line,
                   char *alt_label)
{
    char aname[PDB_ATOM_NAME_STRL + 1], rname[PDB_ATOM_RES_NAME_STRL + 1],
        rnumber[PDB_ATOM_RES_NUMBER_STRL + 1], symbol[PDB_ATOM_SYMBOL_STRL + 1];
    chain_label_t chain_label;
    int flag;
    struct atom *a;

    assert(line);

    if (alt_label) *alt_label = freesasa_pdb_get_alt_coord_label(line);

    freesasa_pdb_get_atom_name(aname, line);
    freesasa_pdb_get_res_name(rname, line);
    freesasa_pdb_get_res_number(rnumber, line);
    chain_label[0] = freesasa_pdb_get_chain_label(line);
    chain_label[1] = '\0';

    flag = freesasa_pdb_get_symbol(symbol, line);
    if (flag == FREESASA_FAIL || (symbol[0] == ' ' && symbol[1] == ' ')) {
        guess_symbol(symbol, aname);
    }

    a = atom_new(rname, rnumber, aname, symbol, chain_label);

    if (a != NULL) {
        a->line = strdup(line);
        if (a->line == NULL) {
            mem_fail();
            atom_free(a);
            a = NULL;
        }
    }

    return a;
}

static struct residues
residues_init(void)
{
    struct residues res;

    res.n = 0;
    res.n_alloc = 0;
    res.first_atom = 0;
    res.reference_area = 0;

    return res;
}

static int
residues_alloc(struct residues *residues)
{
    int new_size;
    void *fa, *ra;

    assert(residues);
    assert(residues->n <= residues->n_alloc);

    if (residues->n == residues->n_alloc) {
        new_size = residues->n_alloc + RESIDUES_CHUNK;
        fa = residues->first_atom;
        ra = residues->reference_area;

        residues->first_atom = realloc(residues->first_atom,
                                       sizeof(int) * new_size);
        if (residues->first_atom == NULL) {
            residues->first_atom = fa;
            return mem_fail();
        }

        residues->reference_area = realloc(residues->reference_area,
                                           sizeof(freesasa_nodearea *) * new_size);
        if (residues->reference_area == NULL) {
            residues->reference_area = ra;
            return mem_fail();
        }

        residues->n_alloc = new_size;
    }
    ++residues->n;
    return FREESASA_SUCCESS;
}

static void
residues_dealloc(struct residues *residues)
{
    int i;

    if (residues) {
        free(residues->first_atom);
        if (residues->reference_area) {
            for (i = 0; i < residues->n; ++i) {
                free(residues->reference_area[i]);
            }
        }
        free(residues->reference_area);
        *residues = residues_init();
    }
}

static struct chains
chains_init(void)
{
    struct chains ch;

    ch.n = 0;
    ch.n_alloc = 0;
    ch.first_atom = NULL;
    ch.labels = NULL;
    ch.short_labels = NULL;

    return ch;
}

static int
chains_alloc(struct chains *chains)
{
    int new_size;
    void *fa, *lbl, *short_lbl;

    assert(chains);
    assert(chains->n <= chains->n_alloc);

    if (chains->n == chains->n_alloc) {
        new_size = chains->n_alloc + CHAINS_CHUNK;
        fa = chains->first_atom;
        lbl = chains->labels;
        short_lbl = chains->short_labels;

        chains->first_atom = realloc(chains->first_atom,
                                     sizeof(int) * new_size);
        if (chains->first_atom == NULL) {
            chains->first_atom = fa;
            return mem_fail();
        }

        chains->labels = realloc(chains->labels, sizeof(chain_label_t) * new_size);
        if (chains->labels == NULL) {
            chains->labels = lbl;
            return mem_fail();
        }

        chains->short_labels = realloc(chains->short_labels, new_size + 1);
        if (chains->short_labels == NULL) {
            chains->short_labels = short_lbl;
            return mem_fail();
        }

        chains->n_alloc = new_size;
    }
    ++chains->n;
    return FREESASA_SUCCESS;
}

static void
chains_dealloc(struct chains *chains)
{
    if (chains) {
        free(chains->first_atom);
        free(chains->labels);
        free(chains->short_labels);
        *chains = chains_init();
    }
}

freesasa_structure *
freesasa_structure_new(void)
{
    freesasa_structure *s = malloc(sizeof(freesasa_structure));

    if (s == NULL) goto memerr;

    s->atoms = atoms_init();
    s->residues = residues_init();
    s->chains = chains_init();
    s->xyz = freesasa_coord_new();
    s->model = 1;
    s->classifier_name = NULL;
    s->cif_ref = 0;
    s->release_cif_ref = 0;

    if (s->xyz == NULL) goto memerr;

    return s;
memerr:
    freesasa_structure_free(s);
    mem_fail();
    return NULL;
}

void freesasa_structure_free(freesasa_structure *s)
{
    if (s != NULL) {
        atoms_dealloc(&s->atoms);
        residues_dealloc(&s->residues);
        chains_dealloc(&s->chains);

        if (s->xyz != NULL) {
            freesasa_coord_free(s->xyz);
        }

        free(s->classifier_name);

        if (s->cif_ref > 0 && s->release_cif_ref != NULL) {
            s->release_cif_ref(s->cif_ref);
        }

        free(s);
    }
}

/**
    This function is called when either the symbol field is missing
    from an ATOM record, or when an atom is added directly using
    freesasa_structure_add_atom() or
    freesasa_structure_add_atom_wopt(). The symbol is in turn only
    used if the atom cannot be recognized by the classifier.
*/
static int
guess_symbol(char *symbol,
             const char *name)
{
    /* if the first position is empty, or a number, assume that it is
       a one letter element e.g. " C ", or "1H " */
    if (name[0] == ' ' || (name[0] >= '1' && name[0] <= '9')) {
        symbol[0] = ' ';
        symbol[1] = name[1];
        symbol[2] = '\0';
    } else {
        /* if the string has padding to the right, it's a
           two-letter element, e.g. "FE  " */
        if (name[3] == ' ') {
            strncpy(symbol, name, 2);
            symbol[2] = '\0';
        } else {
            /* If it's a four-letter string, it's hard to say,
               assume only the first letter signifies the element */
            symbol[0] = ' ';
            symbol[1] = name[0];
            symbol[2] = '\0';
            return freesasa_warn("guessing that atom '%s' is symbol '%s'",
                                 name, symbol);
        }
    }
    return FREESASA_SUCCESS;
}

int structure_has_chain(freesasa_structure *s, const chain_label_t chain_label)
{
    int i;
    for (i = 0; i < s->chains.n; ++i) {
        if (strncmp(s->chains.labels[i], chain_label, sizeof(chain_label_t)) == 0) {
            return 1;
        }
    }
    return 0;
}

static int
structure_add_chain(freesasa_structure *s,
                    const chain_label_t chain_label,
                    int i_latest_atom)
{
    int n;
    if (s->chains.n == 0 || structure_has_chain(s, chain_label) == 0) {

        if (chains_alloc(&s->chains) == FREESASA_FAIL)
            return fail_msg("");

        n = s->chains.n;
        snprintf(s->chains.labels[n - 1], sizeof(chain_label_t), "%s", chain_label);
        s->chains.short_labels[n - 1] = chain_label[0];
        s->chains.short_labels[n] = '\0';

        s->chains.first_atom[n - 1] = i_latest_atom;
    }
    return FREESASA_SUCCESS;
}

static int
structure_add_residue(freesasa_structure *s,
                      const freesasa_classifier *classifier,
                      const struct atom *a,
                      int i_latest_atom)
{
    int n = s->residues.n + 1;
    const freesasa_nodearea *reference = NULL;

    /* register a new residue if it's the first atom, or if the
       residue number or chain label of the current atom is different
       from the previous one */
    if (!(s->residues.n == 0 ||
          (i_latest_atom > 0 &&
           (strcmp(a->res_number, s->atoms.atom[i_latest_atom - 1]->res_number) ||
            strcmp(a->chain_label, s->atoms.atom[i_latest_atom - 1]->chain_label))))) {
        return FREESASA_SUCCESS;
    }

    if (residues_alloc(&s->residues) == FREESASA_FAIL) {
        return fail_msg("");
    }
    s->residues.first_atom[n - 1] = i_latest_atom;

    s->residues.reference_area[n - 1] = NULL;
    reference = freesasa_classifier_residue_reference(classifier, a->res_name);
    if (reference != NULL) {
        s->residues.reference_area[n - 1] = malloc(sizeof(freesasa_nodearea));
        if (s->residues.reference_area[n - 1] == NULL)
            return mem_fail();
        *s->residues.reference_area[n - 1] = *reference;
    }

    return FREESASA_SUCCESS;
}

/**
    Get the radius of an atom, and fail, warn and/or guess depending
    on the options.
 */
static int
structure_check_atom_radius(double *radius,
                            struct atom *a,
                            const freesasa_classifier *classifier,
                            int options)
{
    *radius = freesasa_classifier_radius(classifier, a->res_name, a->atom_name);
    if (*radius < 0) {
        if (options & FREESASA_HALT_AT_UNKNOWN) {
            return fail_msg("atom '%s %s' unknown",
                            a->res_name, a->atom_name);
        } else if (options & FREESASA_SKIP_UNKNOWN) {
            return freesasa_warn("skipping unknown atom '%s %s'",
                                 a->res_name, a->atom_name, a->symbol, *radius);
        } else {
            *radius = freesasa_guess_radius(a->symbol);
            if (*radius < 0) {
                *radius = +0.;
                freesasa_warn("atom '%s %s' unknown and "
                              "can't guess radius of symbol '%s', "
                              "assigning radius 0 A",
                              a->res_name, a->atom_name, a->symbol);
            } else {
                freesasa_warn("atom '%s %s' unknown, guessing element is '%s', "
                              "and radius %.3f A",
                              a->res_name, a->atom_name, a->symbol, *radius);
            }
            /* do not return FREESASA_WARN here, because we will keep the atom */
        }
    }
    return FREESASA_SUCCESS;
}

static int
structure_register_classifier(freesasa_structure *structure,
                              const freesasa_classifier *classifier)
{
    if (structure->classifier_name == NULL) {
        structure->classifier_name = strdup(freesasa_classifier_name(classifier));
        if (structure->classifier_name == NULL) {
            return mem_fail();
        }
    } else if (strcmp(structure->classifier_name, freesasa_classifier_name(classifier)) != 0) {
        structure->classifier_name = strdup(FREESASA_CONFLICTING_CLASSIFIERS);
        if (structure->classifier_name == NULL) {
            return mem_fail();
        }
        return FREESASA_WARN;
    }

    return FREESASA_SUCCESS;
}

/**
   Adds an atom to the structure using the rules specified by
   'options'. If it includes FREESASA_RADIUS_FROM_* a dummy radius is
   assigned and the caller is expected to replace it with a correct
   radius later.

   The atom a should be a pointer to a heap address, this will not be cloned.
 */
static int
structure_add_atom(freesasa_structure *structure,
                   struct atom *atom,
                   double *xyz,
                   const freesasa_classifier *classifier,
                   int options)
{
    int na, ret;
    double r;

    assert(structure);
    assert(atom);
    assert(xyz);

    /* let the stricter option override if both are specified */
    if (options & FREESASA_SKIP_UNKNOWN && options & FREESASA_HALT_AT_UNKNOWN)
        options &= ~FREESASA_SKIP_UNKNOWN;

    if (classifier == NULL) {
        classifier = &freesasa_default_classifier;
    }
    structure_register_classifier(structure, classifier);

    /* calculate radius and check if we should keep the atom (based on options) */
    if (options & FREESASA_RADIUS_FROM_OCCUPANCY) {
        r = 1; /* fix it later */
    } else {
        ret = structure_check_atom_radius(&r, atom, classifier, options);
        if (ret == FREESASA_FAIL) return fail_msg("halting at unknown atom");
        if (ret == FREESASA_WARN) return FREESASA_WARN;
    }
    assert(r >= 0);

    /* If it's a keeper, allocate memory */
    if (atoms_alloc(&structure->atoms) == FREESASA_FAIL)
        return fail_msg("");
    na = structure->atoms.n;

    /* Store coordinates */
    if (freesasa_coord_append(structure->xyz, xyz, 1) == FREESASA_FAIL)
        return mem_fail();

    /* Check if this is a new chain and if so add it */
    if (structure_add_chain(structure, atom->chain_label, na - 1) == FREESASA_FAIL)
        return mem_fail();

    /* Check if this is a new residue, and if so add it */
    if (structure_add_residue(structure, classifier, atom, na - 1) == FREESASA_FAIL)
        return mem_fail();

    atom->the_class = freesasa_classifier_class(classifier, atom->res_name, atom->atom_name);
    atom->res_index = structure->residues.n - 1;
    structure->atoms.radius[na - 1] = r;
    structure->atoms.atom[na - 1] = atom;

    return FREESASA_SUCCESS;
}

/**
    Handles the reading of PDB-files, returns NULL if problems reading
    or input or malloc failure. Error-messages should explain what
    went wrong.
 */
static freesasa_structure *
from_pdb_impl(FILE *pdb_file,
              struct file_range it,
              const freesasa_classifier *classifier,
              int options)
{
    char line[PDB_MAX_LINE_STRL];
    char alt, the_alt = ' ';
    double v[3], r;
    int ret;
    struct atom *a = NULL;
    freesasa_structure *s = freesasa_structure_new();

    assert(pdb_file);

    if (s == NULL) return NULL;

    fseek(pdb_file, it.begin, SEEK_SET);

    while (fgets(line, PDB_MAX_LINE_STRL, pdb_file) != NULL && ftell(pdb_file) <= it.end) {

        if (strncmp("ATOM", line, 4) == 0 || ((options & FREESASA_INCLUDE_HETATM) &&
                                              (strncmp("HETATM", line, 6) == 0))) {
            if (freesasa_pdb_ishydrogen(line) &&
                !(options & FREESASA_INCLUDE_HYDROGEN))
                continue;

            a = atom_new_from_line(line, &alt);
            if (a == NULL)
                goto cleanup;

            if ((alt != ' ' && the_alt == ' ') || (alt == ' '))
                the_alt = alt;
            else if (alt != ' ' && alt != the_alt) {
                atom_free(a);
                a = NULL;
                continue;
            }

            ret = freesasa_pdb_get_coord(v, line);
            if (ret == FREESASA_FAIL)
                goto cleanup;

            ret = structure_add_atom(s, a, v, classifier, options);
            if (ret == FREESASA_FAIL) {
                goto cleanup;
            } else if (ret == FREESASA_WARN) {
                atom_free(a);
                a = NULL;
                continue;
            }

            if (options & FREESASA_RADIUS_FROM_OCCUPANCY) {
                ret = freesasa_pdb_get_occupancy(&r, line);
                if (ret == FREESASA_FAIL)
                    goto cleanup;
                s->atoms.radius[s->atoms.n - 1] = r;
            }
        }

        if (!(options & FREESASA_JOIN_MODELS)) {
            if (strncmp("MODEL", line, 5) == 0) sscanf(line + 10, "%d", &s->model);
            if (strncmp("ENDMDL", line, 6) == 0) break;
        }
    }

    if (s->atoms.n == 0) {
        fail_msg("input had no valid ATOM or HETATM lines");
        goto cleanup;
    }

    return s;

cleanup:
    fail_msg("");
    atom_free(a);
    freesasa_structure_free(s);
    return NULL;
}

static int
structure_add_atom_wopt_impl(freesasa_structure *structure,
                             const char *atom_name,
                             const char *residue_name,
                             const char *residue_number,
                             const char *symbol,
                             const char *chain_label,
                             double x, double y, double z,
                             const freesasa_classifier *classifier,
                             int options)
{
    struct atom *a;
    char my_symbol[PDB_ATOM_SYMBOL_STRL + 1];
    double v[3] = {x, y, z};
    int ret, warn = 0;

    assert(structure);
    assert(atom_name);
    assert(residue_name);
    assert(residue_number);
    assert(chain_label);

    /* this option can not be used here, and needs to be unset */
    options &= ~FREESASA_RADIUS_FROM_OCCUPANCY;

    if (symbol != NULL) {
        strncpy(my_symbol, symbol, sizeof(my_symbol));
    } else if (guess_symbol(my_symbol, atom_name) == FREESASA_WARN &&
               options & FREESASA_SKIP_UNKNOWN) {
        ++warn;
    }

    a = atom_new(residue_name, residue_number, atom_name, my_symbol, chain_label);
    if (a == NULL) return mem_fail();

    ret = structure_add_atom(structure, a, v, classifier, options);

    if (ret == FREESASA_FAIL ||
        (ret == FREESASA_WARN && options & FREESASA_SKIP_UNKNOWN))
        atom_free(a);

    if (!ret && warn) return FREESASA_WARN;

    return ret;
}

int freesasa_structure_add_atom_wopt(freesasa_structure *structure,
                                     const char *atom_name,
                                     const char *residue_name,
                                     const char *residue_number,
                                     char chain_label,
                                     double x, double y, double z,
                                     const freesasa_classifier *classifier,
                                     int options)
{
    chain_label_t my_chain_label = {chain_label, '\0'};

    return structure_add_atom_wopt_impl(structure, atom_name, residue_name, residue_number, NULL,
                                        my_chain_label, x, y, z, classifier, options);
}

int freesasa_structure_add_atom(freesasa_structure *structure,
                                const char *atom_name,
                                const char *residue_name,
                                const char *residue_number,
                                char chain_label,
                                double x, double y, double z)
{
    chain_label_t my_chain_label = {chain_label, '\0'};

    return structure_add_atom_wopt_impl(structure, atom_name, residue_name, residue_number, NULL,
                                        my_chain_label, x, y, z, NULL, 0);
}

int freesasa_structure_add_cif_atom(freesasa_structure *structure,
                                    freesasa_cif_atom *atom,
                                    const freesasa_classifier *classifier,
                                    int options)
{
    char res_number[PDB_ATOM_RES_NUMBER_STRL + 1];

    if (atom->pdbx_PDB_ins_code[0] != '?') {
        snprintf(res_number, sizeof res_number, "%s%c", atom->auth_seq_id, atom->pdbx_PDB_ins_code[0]);
    } else {
        snprintf(res_number, sizeof res_number, "%s", atom->auth_seq_id);
    }

    chain_label_t chain_label = {atom->auth_asym_id, '\0'};

    return structure_add_atom_wopt_impl(structure, atom->auth_atom_id, atom->auth_comp_id,
                                        res_number, atom->type_symbol, chain_label,
                                        atom->Cartn_x, atom->Cartn_y, atom->Cartn_z,
                                        classifier, options);
}

int freesasa_structure_add_cif_atom_lcl(freesasa_structure *structure,
                                        freesasa_cif_atom_lcl *atom,
                                        const freesasa_classifier *classifier,
                                        int options)
{
    char res_number[PDB_ATOM_RES_NUMBER_STRL + 1];

    if (atom->pdbx_PDB_ins_code[0] != '?') {
        snprintf(res_number, sizeof res_number, "%s%c", atom->auth_seq_id, atom->pdbx_PDB_ins_code[0]);
    } else {
        snprintf(res_number, sizeof res_number, "%s", atom->auth_seq_id);
    }

    return structure_add_atom_wopt_impl(structure, atom->auth_atom_id, atom->auth_comp_id,
                                        res_number, atom->type_symbol, atom->auth_asym_id,
                                        atom->Cartn_x, atom->Cartn_y, atom->Cartn_z,
                                        classifier, options);
}

freesasa_structure *
freesasa_structure_from_pdb(FILE *pdb_file,
                            const freesasa_classifier *classifier,
                            int options)
{
    assert(pdb_file);
    return from_pdb_impl(pdb_file, freesasa_whole_file(pdb_file),
                         classifier, options);
}

freesasa_structure **
freesasa_structure_array(FILE *pdb,
                         int *n,
                         const freesasa_classifier *classifier,
                         int options)
{
    struct file_range *models = NULL, *chains = NULL;
    struct file_range whole_file;
    int n_models = 0, n_chains = 0, j0, n_new_chains, i, j;
    freesasa_structure **ss = NULL, **ssb;

    assert(pdb);
    assert(n);

    if (!(options & FREESASA_SEPARATE_MODELS ||
          options & FREESASA_SEPARATE_CHAINS)) {
        fail_msg("options need to specify at least one of FREESASA_SEPARATE_CHAINS "
                 "and FREESASA_SEPARATE_MODELS");
        return NULL;
    }

    whole_file = freesasa_whole_file(pdb);
    n_models = freesasa_pdb_get_models(pdb, &models);

    if (n_models == FREESASA_FAIL) {
        fail_msg("problems reading PDB-file");
        return NULL;
    }
    if (n_models == 0) {
        models = &whole_file;
        n_models = 1;
    }

    /* only keep first model if option not provided */
    if (!(options & FREESASA_SEPARATE_MODELS)) n_models = 1;

    /* for each model read chains if requested */
    if (options & FREESASA_SEPARATE_CHAINS) {
        for (i = 0; i < n_models; ++i) {
            chains = NULL;
            n_new_chains = freesasa_pdb_get_chains(pdb, models[i], &chains, options);

            if (n_new_chains == FREESASA_FAIL) goto cleanup;
            if (n_new_chains == 0) {
                freesasa_warn("in %s(): no chains found (in model %d)", __func__, i + 1);
                continue;
            }

            ssb = ss;
            ss = realloc(ss, sizeof(freesasa_structure *) * (n_chains + n_new_chains));

            if (!ss) {
                ss = ssb;
                mem_fail();
                goto cleanup;
            }

            j0 = n_chains;
            n_chains += n_new_chains;

            for (j = 0; j < n_new_chains; ++j)
                ss[j0 + j] = NULL;

            for (j = 0; j < n_new_chains; ++j) {
                ss[j0 + j] = from_pdb_impl(pdb, chains[j], classifier, options);
                if (ss[j0 + j] == NULL) goto cleanup;
                ss[j0 + j]->model = i + 1;
            }

            free(chains);
        }
        *n = n_chains;
    } else {
        ss = malloc(sizeof(freesasa_structure *) * n_models);
        if (!ss) {
            mem_fail();
            goto cleanup;
        }

        for (i = 0; i < n_models; ++i)
            ss[i] = NULL;
        *n = n_models;

        for (i = 0; i < n_models; ++i) {
            ss[i] = from_pdb_impl(pdb, models[i], classifier, options);
            if (ss[i] == NULL) goto cleanup;
            ss[i]->model = i + 1;
        }
    }

    if (*n == 0) goto cleanup;

    if (models != &whole_file) free(models);

    return ss;

cleanup:
    if (ss)
        for (i = 0; i < *n; ++i)
            freesasa_structure_free(ss[i]);
    if (models != &whole_file) free(models);
    free(chains);
    *n = 0;
    free(ss);
    return NULL;
}

freesasa_structure *
freesasa_structure_get_chains(const freesasa_structure *structure,
                              const char *chains,
                              const freesasa_classifier *classifier,
                              int options)
{
    freesasa_structure *new_s;
    struct atom *ai;
    int i, res;
    char *c;
    const double *v;

    assert(structure);
    if (strlen(chains) == 0) return NULL;

    new_s = freesasa_structure_new();

    if (!new_s) {
        mem_fail();
        return NULL;
    }

    new_s->model = structure->model;

    for (i = 0; i < structure->atoms.n; ++i) {
        ai = structure->atoms.atom[i];
        c = ai->chain_label;
        if (strchr(chains, c[0]) != NULL) {
            v = freesasa_coord_i(structure->xyz, i);
            res = structure_add_atom_wopt_impl(new_s, ai->atom_name,
                                               ai->res_name, ai->res_number, ai->symbol,
                                               c, v[0], v[1], v[2], classifier, options);
            if (res == FREESASA_FAIL) {
                fail_msg("");
                goto cleanup;
            }
        }
    }

    /* the following two tests could have been done by comparing the
       chain-strings before the loop, but this logic is simpler. */
    if (new_s->atoms.n == 0) {
        goto cleanup;
    }
    if (new_s->chains.n != strlen(chains)) {
        fail_msg("structure has chains '%s', but '%s' requested",
                 structure->chains.labels, chains);
        goto cleanup;
    }

    return new_s;

cleanup:
    freesasa_structure_free(new_s);
    return NULL;
}

static int chain_group_has_chain(const freesasa_chain_group *chains, chain_label_t chain)
{
    assert(chains);

    for (size_t i = 0; i < chains->n; ++i) {
        if (strncmp(chains->chains[i], chain, sizeof(chain_label_t)) == 0) {
            return 1;
        }
    }

    return 0;
}

freesasa_structure *
freesasa_structure_get_chains_lcl(const freesasa_structure *structure,
                                  const freesasa_chain_group *chains,
                                  const freesasa_classifier *classifier,
                                  int options)
{
    freesasa_structure *new_s;
    struct atom *ai;
    int i, res;
    char *c;
    const double *v;

    assert(structure);
    assert(chains);
    if (chains->n == 0) return NULL;

    new_s = freesasa_structure_new();

    if (!new_s) {
        mem_fail();
        return NULL;
    }

    new_s->model = structure->model;

    for (i = 0; i < structure->atoms.n; ++i) {
        ai = structure->atoms.atom[i];
        c = ai->chain_label;
        if (chain_group_has_chain(chains, c)) {
            v = freesasa_coord_i(structure->xyz, i);
            res = structure_add_atom_wopt_impl(new_s, ai->atom_name,
                                               ai->res_name, ai->res_number, ai->symbol,
                                               c, v[0], v[1], v[2], classifier, options);
            if (res == FREESASA_FAIL) {
                fail_msg("");
                goto cleanup;
            }
        }
    }

    /* the following two tests could have been done by comparing the
       chain-strings before the loop, but this logic is simpler. */
    if (new_s->atoms.n == 0) {
        goto cleanup;
    }
    if (new_s->chains.n != chains->n) {
        fail_msg("structure has chains '%s', but '%s' requested",
                 structure->chains.labels, chains);
        goto cleanup;
    }

    return new_s;

cleanup:
    freesasa_structure_free(new_s);
    return NULL;
}

const char *
freesasa_structure_chain_labels(const freesasa_structure *structure)
{
    assert(structure);
    return structure->chains.short_labels;
}

int freesasa_structure_number_chains(const freesasa_structure *structure)
{
    assert(structure);

    return structure->chains.n;
}

const char *
freesasa_structure_chain_label(const freesasa_structure *structure, int index)
{
    assert(structure);
    assert(index >= 0 && index < structure->chains.n);

    return structure->chains.labels[index];
}

const coord_t *
freesasa_structure_xyz(const freesasa_structure *structure)
{
    assert(structure);
    return structure->xyz;
}

int freesasa_structure_n(const freesasa_structure *structure)
{
    assert(structure);
    return structure->atoms.n;
}

int freesasa_structure_n_residues(const freesasa_structure *structure)
{
    assert(structure);
    return structure->residues.n;
}

const char *
freesasa_structure_atom_name(const freesasa_structure *structure,
                             int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->atom_name;
}

const char *
freesasa_structure_atom_res_name(const freesasa_structure *structure,
                                 int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->res_name;
}

const char *
freesasa_structure_atom_res_number(const freesasa_structure *structure,
                                   int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->res_number;
}

char freesasa_structure_atom_chain(const freesasa_structure *structure,
                                   int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->chain_label[0];
}

const char *
freesasa_structure_atom_chain_lcl(const freesasa_structure *structure,
                                  int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->chain_label;
}

const char *
freesasa_structure_atom_symbol(const freesasa_structure *structure,
                               int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->symbol;
}

double
freesasa_structure_atom_radius(const freesasa_structure *structure,
                               int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.radius[i];
}

void freesasa_structure_atom_set_radius(freesasa_structure *structure,
                                        int i,
                                        double radius)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    structure->atoms.radius[i] = radius;
}

freesasa_atom_class
freesasa_structure_atom_class(const freesasa_structure *structure,
                              int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->the_class;
}

const char *
freesasa_structure_atom_pdb_line(const freesasa_structure *structure,
                                 int i)
{
    assert(structure);
    assert(i < structure->atoms.n && i >= 0);
    return structure->atoms.atom[i]->line;
}
const freesasa_nodearea *
freesasa_structure_residue_reference(const freesasa_structure *structure,
                                     int r_i)
{
    assert(structure);
    assert(r_i >= 0 && r_i < structure->residues.n);

    return structure->residues.reference_area[r_i];
}
int freesasa_structure_residue_atoms(const freesasa_structure *structure,
                                     int r_i,
                                     int *first,
                                     int *last)
{
    int naa;
    assert(structure);
    assert(first);
    assert(last);

    naa = structure->residues.n;

    assert(r_i >= 0 && r_i < naa);

    *first = structure->residues.first_atom[r_i];
    if (r_i == naa - 1)
        *last = structure->atoms.n - 1;
    else
        *last = structure->residues.first_atom[r_i + 1] - 1;
    assert(*last >= *first);

    return FREESASA_SUCCESS;
}

const char *
freesasa_structure_residue_name(const freesasa_structure *structure,
                                int r_i)
{
    assert(structure);
    assert(r_i < structure->residues.n && r_i >= 0);
    return structure->atoms.atom[structure->residues.first_atom[r_i]]->res_name;
}

const char *
freesasa_structure_residue_number(const freesasa_structure *structure,
                                  int r_i)
{
    assert(structure);
    assert(r_i < structure->residues.n && r_i >= 0);
    return structure->atoms.atom[structure->residues.first_atom[r_i]]->res_number;
}

char freesasa_structure_residue_chain(const freesasa_structure *structure,
                                      int r_i)
{
    assert(structure);
    assert(r_i < structure->residues.n && r_i >= 0);

    return structure->atoms.atom[structure->residues.first_atom[r_i]]->chain_label[0];
}

const char *
freesasa_structure_residue_chain_lcl(const freesasa_structure *structure,
                                     int r_i)
{
    assert(structure);
    assert(r_i < structure->residues.n && r_i >= 0);

    return structure->atoms.atom[structure->residues.first_atom[r_i]]->chain_label;
}

int freesasa_structure_n_chains(const freesasa_structure *structure)
{
    return structure->chains.n;
}

int freesasa_structure_chain_index_lcl(const freesasa_structure *structure,
                                       const char *chain)
{
    int i;

    assert(structure);

    for (i = 0; i < structure->chains.n; ++i) {
        if (strncmp(structure->chains.labels[i], chain, sizeof(chain_label_t)) == 0) {
            return i;
        }
    }

    return fail_msg("chain '%s' not found", chain);
}

/** Not public */
int freesasa_structure_chain_index(const freesasa_structure *structure,
                                   char chain)
{
    chain_label_t chain_label = {chain, '\0'};
    return freesasa_structure_chain_index_lcl(structure, chain_label);
}

int freesasa_structure_chain_atoms_lcl(const freesasa_structure *structure,
                                       const char *chain,
                                       int *first,
                                       int *last)
{
    int c_i, n;

    assert(structure);

    c_i = freesasa_structure_chain_index_lcl(structure, chain);
    n = freesasa_structure_n_chains(structure);

    if (c_i < 0) return fail_msg("");

    *first = structure->chains.first_atom[c_i];
    if (c_i == n - 1)
        *last = freesasa_structure_n(structure) - 1;
    else
        *last = structure->chains.first_atom[c_i + 1] - 1;
    assert(*last >= *first);

    return FREESASA_SUCCESS;
}

int freesasa_structure_chain_atoms(const freesasa_structure *structure,
                                   char chain,
                                   int *first,
                                   int *last)
{
    chain_label_t chain_label = {chain, '\0'};
    return freesasa_structure_chain_atoms_lcl(structure, chain_label, first, last);
}

int freesasa_structure_chain_residues_lcl(const freesasa_structure *structure,
                                          const char *chain,
                                          int *first,
                                          int *last)
{
    int first_atom, last_atom;

    assert(structure);

    if (freesasa_structure_chain_atoms_lcl(structure, chain, &first_atom, &last_atom))
        return fail_msg("");

    *first = structure->atoms.atom[first_atom]->res_index;
    *last = structure->atoms.atom[last_atom]->res_index;

    return FREESASA_SUCCESS;
}

int freesasa_structure_chain_residues(const freesasa_structure *structure,
                                      char chain,
                                      int *first,
                                      int *last)
{
    chain_label_t chain_label = {chain, '\0'};
    return freesasa_structure_chain_residues_lcl(structure, chain_label, first, last);
}

const char *
freesasa_structure_classifier_name(const freesasa_structure *structure)
{
    assert(structure);
    return structure->classifier_name;
}

int freesasa_structure_model(const freesasa_structure *structure)
{
    return structure->model;
}

const double *
freesasa_structure_coord_array(const freesasa_structure *structure)
{
    return freesasa_coord_all(structure->xyz);
}

const double *
freesasa_structure_radius(const freesasa_structure *structure)
{
    assert(structure);
    return structure->atoms.radius;
}

void freesasa_structure_set_radius(freesasa_structure *structure,
                                   const double *radii)
{
    assert(structure);
    assert(radii);
    memcpy(structure->atoms.radius, radii, structure->atoms.n * sizeof(double));
}

void freesasa_structure_set_model(freesasa_structure *structure,
                                  int model)
{
    structure->model = model;
}

void freesasa_structure_set_cif_ref(freesasa_structure *structure,
                                    size_t ref,
                                    void (*release_func)(size_t))
{
    assert(structure);
    structure->cif_ref = ref;
    structure->release_cif_ref = release_func;
}

size_t freesasa_structure_cif_ref(const freesasa_structure *structure)
{
    assert(structure);
    return structure->cif_ref;
}
