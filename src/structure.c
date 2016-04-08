#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "freesasa_internal.h"
#include "pdb.h"
#include "classifier.h"
#include "coord.h"

struct atom {
    char *res_name;
    char *res_number;
    char *atom_name;
    char *symbol;
    char *descriptor;
    char *line;
    char chain_label;
};

static const struct atom empty_atom =
     {NULL, NULL, NULL, NULL, NULL, NULL, '\0'};

struct freesasa_structure {
    struct atom **a;
    coord_t *xyz;
    double *radius;
    int number_atoms;
    int number_residues;
    int number_chains;
    int model; // model number
    char *chains; //all chain labels found (as string)
    int *res_first_atom; // first atom of each residue
    int *chain_first_atom; // first atom of each chain
    char **res_desc;
};

static const struct freesasa_structure empty_structure = 
    {NULL,NULL,NULL,0,0,0,0,NULL,NULL,NULL,NULL};

static int
guess_symbol(char *symbol,
             const char *name);

static void
atom_free(struct atom *a)
{
    if (a == NULL) return;
    free(a->res_name);
    free(a->res_number);
    free(a->atom_name);
    free(a->symbol);
    free(a->descriptor);
    free(a->line);
    free(a);
}

static struct atom *
atom_new(const char *residue_name,
         const char *residue_number,
         const char *atom_name,
         const char *symbol,
         char chain_label)
{
    struct atom *a = malloc(sizeof(struct atom));
    int dlen = 0;
    if (a == NULL) { 
        mem_fail(); 
        return NULL; 
    }

    *a = empty_atom;

    a->line = NULL;
    a->chain_label = chain_label;

    a->res_name = strdup(residue_name);
    a->res_number = strdup(residue_number);
    a->atom_name = strdup(atom_name);
    a->symbol = strdup(symbol);

    dlen = strlen(residue_number) + strlen(residue_name)
        + strlen(atom_name) + 4;
    a->descriptor = malloc(dlen+1);

    if (!a->res_name || !a->res_number || !a->atom_name ||
        !a->symbol || !a->descriptor) {
        mem_fail();
        atom_free(a);
        return NULL;
    }

    sprintf(a->descriptor,"%c %s %s %s",
            chain_label,residue_number,residue_name,atom_name);

    return a;
}

static struct atom *
atom_new_from_line(const char *line,
                   char *alt_label) 
{
    assert(line);
    const int buflen = strlen(line);
    int flag;
    struct atom *a;
    char aname[buflen], rname[buflen], rnumber[buflen], symbol[buflen];

    if (alt_label) *alt_label = freesasa_pdb_get_alt_coord_label(line);

    freesasa_pdb_get_atom_name(aname, line);
    freesasa_pdb_get_res_name(rname, line);
    freesasa_pdb_get_res_number(rnumber, line);

    flag = freesasa_pdb_get_symbol(symbol, line);
    if (flag == FREESASA_FAIL) guess_symbol(symbol,aname);

    a = atom_new(rname, rnumber, aname, symbol, freesasa_pdb_get_chain_label(line));
    
    if (a == NULL) return NULL;

    a->line = strdup(line);
    
    if (a->line == NULL) {
        mem_fail();
        atom_free(a);
        return NULL;
    }

    return a;
}

freesasa_structure*
freesasa_structure_new(void)
{
    freesasa_structure *s = malloc(sizeof(freesasa_structure));

    if (s == NULL) {
        mem_fail();
        return NULL;
    }

    *s = empty_structure;

    s->chains = malloc(1);
    s->a = malloc(sizeof(struct atom*));
    s->xyz = freesasa_coord_new();

    if (!s->chains || !s->a || !s->xyz) {
        freesasa_structure_free(s);
        mem_fail();
        return NULL;
    }

    s->chains[0] = '\0';

    return s;
}

void
freesasa_structure_free(freesasa_structure *s)
{
    if (s == NULL) return;
    if (s->a) {
        for (int i = 0; i < s->number_atoms; ++i) 
            if (s->a[i]) atom_free(s->a[i]);
        free(s->a);
    }
    if (s->xyz) freesasa_coord_free(s->xyz);
    if (s->res_desc) {
        for (int i = 0; i < s->number_residues; ++i)
            if (s->res_desc[i]) free(s->res_desc[i]);
       free(s->res_desc);
    }
    free(s->radius);
    free(s->res_first_atom);
    free(s->chain_first_atom);
    free(s->chains);
    free(s);
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
    // if the first position is empty, assume that it is a one letter element
    // e.g. " C  "
    if (name[0] == ' ') { 
        symbol[0] = ' ';
        symbol[1] = name[1];
        symbol[2] = '\0';
    } else { 
        // if the string has padding to the right, it's a
        // two-letter element, e.g. "FE  "
        if (name[3] == ' ') {
            strncpy(symbol,name,2);
            symbol[2] = '\0';
        } else { 
            // If it's a four-letter string, it's hard to say,
            // assume only the first letter signifies the element
            symbol[0] = ' ';
            symbol[1] = name[0];
            symbol[2] = '\0';
            return freesasa_warn("Guessing that atom '%s' is symbol '%s'",
                                 name,symbol);
        }
    }
    return FREESASA_SUCCESS;
}
static int
structure_add_chain(freesasa_structure *s,
                    char chain_label,
                    int i_latest_atom)
{
    if (strchr(s->chains,chain_label) == NULL) {
        char *sc = s->chains;
        int *cfa = s->chain_first_atom;
        int n = ++s->number_chains;
        s->chains = realloc(s->chains,n + 1);
        s->chain_first_atom = realloc(s->chain_first_atom, n*sizeof(int));
        if (s->chains && s->chain_first_atom) {
            s->chains[n-1] = chain_label;
            s->chains[n] = '\0';
            assert (strlen(s->chains) == s->number_chains);
            s->chain_first_atom[n-1] = i_latest_atom;
        } else {
            s->chains = sc;
            s->chain_first_atom = cfa;
            return mem_fail();
        }
    }
    return FREESASA_SUCCESS;
}

static int
structure_add_residue(freesasa_structure *s, const struct atom *a, int i)
{
    int n = s->number_residues+1;
    int *rfa = s->res_first_atom;
    char **rd = s->res_desc;

    rfa = realloc(rfa, sizeof(int) * n);
    if (!rfa) return mem_fail();
    rfa[n-1] = i;
    s->res_first_atom = rfa;

    rd = realloc(rd, sizeof(char*)*n);
    if (!rd) return mem_fail();

    rd[n-1] = malloc(strlen(a->res_number) + strlen(a->res_name)+4);
    if (!rd[n-1]) {
        free(rd);
        return mem_fail();
    }
    sprintf(rd[n-1], "%c %s %s", a->chain_label, a->res_number, a->res_name);

    s->res_desc = rd;
    ++s->number_residues;

    return FREESASA_SUCCESS;
}

/**
    Get the radius of an atom, and fail, warn and/or guess depending
    on the options.
 */
static int
structure_check_atom_radius(double *radius,
                            struct atom *a,
                            const freesasa_classifier* classifier,
                            int options)
{
    *radius = classifier->radius(a->res_name, a->atom_name, classifier);
    if (*radius < 0) {
        if (options & FREESASA_HALT_AT_UNKNOWN) {
            return freesasa_fail("in %s(): atom '%s %s' unknown.",
                                 __func__, a->res_name, a->atom_name);
        } else if (options & FREESASA_SKIP_UNKNOWN) {
            return freesasa_warn("Skipping unknown atom '%s %s'.",
                                 a->res_name, a->atom_name, a->symbol, *radius);
        } else {
            *radius = freesasa_guess_radius(a->symbol);
            if (*radius < 0) {
                *radius = +0.;
                freesasa_warn("Atom '%s %s' unknown and ",
                              "can't guess radius of symbol '%s'. "
                              "Assigning radius 0 A.",
                              a->res_name, a->atom_name, a->symbol);
            } else {
                freesasa_warn("Atom '%s %s' unknown, guessing element is '%s', "
                              "and radius %.3f A.",
                              a->res_name, a->atom_name, a->symbol, *radius);
            }
            // do not return FREESASA_WARN here, because we will keep the atom
        }
    }
    return FREESASA_SUCCESS;
}

/**
   Adds an atom to the structure using the rules specified by
   'options'. If it includes FREESASA_RADIUS_FROM_* a dummy radius is
   assigned and the caller is expected to replace it with a correct
   radius later.
 */
static int
structure_add_atom(freesasa_structure *s,
                   struct atom *a,
                   double *xyz,
                   const freesasa_classifier* classifier,
                   int options)
{
    assert(s); assert(a); assert(xyz);
    int na, ret;
    double r, *pr = s->radius;
    struct atom **pa = s->a;
    
    if (classifier == NULL) {
        classifier = &freesasa_default_classifier;
    }

    // calculate radius and check if we should keep the atom (based on options)
    if (options & FREESASA_RADIUS_FROM_OCCUPANCY) {
        r = 1; // fix it later
    } else {
        ret = structure_check_atom_radius(&r, a, classifier, options);
        if (ret == FREESASA_FAIL) return fail_msg("Halting at unknown atom.");
        if (ret == FREESASA_WARN) return FREESASA_WARN;
    }
    assert(r >= 0);

    // if it's a keeper store the radius
    na = s->number_atoms+1;
    if ((s->radius = realloc(s->radius,sizeof(double)*na)) == NULL) {
        s->radius = pr;
        return mem_fail();
    }
    s->radius[na-1] = r;

    // allocate memory for atom, add chain
    if ((s->a = realloc(s->a,sizeof(struct atom*)*na)) == NULL) {
        s->a = pa;
        return mem_fail();
    }

    if (freesasa_coord_append(s->xyz, xyz, 1)) return mem_fail();
    if (structure_add_chain(s, a->chain_label, na - 1)) return mem_fail();

    /* register a new residue if it's the first atom, or if the
       residue number or chain label of the current atom is different
       from the previous one */
    if ( s->number_residues == 0 ||
         (na > 1 && 
          (strcmp(a->res_number, s->a[na-2]->res_number) ||
           a->chain_label != s->a[na-2]->chain_label) )) {
        if (structure_add_residue(s, a, na - 1)) return mem_fail();
    }

    // by doing this last, we can free as much memory as possible if anything fails
    s->a[na-1] = a;
    ++s->number_atoms;

    return FREESASA_SUCCESS;
}

/**
    Handles the reading of PDB-files, returns NULL if problems reading
    or input or malloc failure. Error-messages should explain what
    went wrong.
 */
static freesasa_structure*
from_pdb_impl(FILE *pdb_file,
              struct file_range it,
              const freesasa_classifier *classifier,
              int options)
{
    assert(pdb_file);
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    char alt, the_alt = ' ';
    double v[3], r;
    struct atom *a = NULL;
    freesasa_structure *s = freesasa_structure_new();
 
    if (s == NULL) return NULL;
    
    fseek(pdb_file,it.begin,SEEK_SET);
    
    while (getline(&line, &len, pdb_file) != -1 && ftell(pdb_file) <= it.end) {
        
        if (strncmp("ATOM",line,4)==0 || ( (options & FREESASA_INCLUDE_HETATM) &&
                                           (strncmp("HETATM", line, 6) == 0) )) {

            if (freesasa_pdb_ishydrogen(line) &&
                !(options & FREESASA_INCLUDE_HYDROGEN))
                continue;

            if (!(a = atom_new_from_line(line, &alt))) goto cleanup;

            if ((alt != ' ' && the_alt == ' ') || (alt == ' '))
                the_alt = alt;
            else if (alt != ' ' && alt != the_alt)
                continue;

            if (freesasa_pdb_get_coord(v, line) == FREESASA_FAIL ||
                structure_add_atom(s, a, v, classifier, options) == FREESASA_FAIL) 
                goto cleanup;

            if (options & FREESASA_RADIUS_FROM_OCCUPANCY) {
                if (freesasa_pdb_get_occupancy(&r, line) == FREESASA_SUCCESS) {
                    s->radius[s->number_atoms-1] = r;
                }
            }
        }

        if (! (options & FREESASA_JOIN_MODELS)) {
            if (strncmp("MODEL",line,5)==0)  sscanf(line+10, "%d", &s->model);
            if (strncmp("ENDMDL",line,6)==0) break;
        }
    }
    
    if (s->number_atoms == 0) {
        freesasa_fail("Input had no valid ATOM or HETATM lines.");
        goto cleanup;
    }

    free(line);
    return s;

 cleanup:
    fail_msg("");
    free(line);
    atom_free(a);
    freesasa_structure_free(s);
    return NULL;
}


int
freesasa_structure_add_atom_wopt(freesasa_structure *s,
                                 const char *atom_name,
                                 const char *residue_name,
                                 const char *residue_number,
                                 char chain_label,
                                 double x, double y, double z,
                                 const freesasa_classifier *classifier,
                                 int options)
{
    assert(s);
    assert(atom_name); assert(residue_name); assert(residue_number);

    struct atom *a;
    char symbol[PDB_ATOM_SYMBOL_STRL+1];
    double v[3] = {x,y,z};
    int ret, warn = 0;

    // this option can not be used here, and needs to be unset
    options &= ~FREESASA_RADIUS_FROM_OCCUPANCY;

    if (guess_symbol(symbol, atom_name) == FREESASA_WARN &&
        options & FREESASA_SKIP_UNKNOWN) 
        ++warn;

    a = atom_new(residue_name, residue_number, atom_name, symbol, chain_label);
    if (a == NULL) return mem_fail();

    ret = structure_add_atom(s, a, v, classifier, options);

    if (ret == FREESASA_FAIL) {
        atom_free(a);
        return ret;
    }

    if (warn) return FREESASA_WARN;

    return ret;
}

int 
freesasa_structure_add_atom(freesasa_structure *s,
                            const char *atom_name,
                            const char *residue_name,
                            const char *residue_number,
                            char chain_label,
                            double x, double y, double z)
{
    return freesasa_structure_add_atom_wopt(s, atom_name, residue_name, residue_number,
                                            chain_label, x, y, z, NULL, 0);
}

freesasa_structure *
freesasa_structure_from_pdb(FILE *pdb_file,
                            const freesasa_classifier* classifier,
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
    assert(pdb);
    assert(n);

    struct file_range *models = NULL, *chains = NULL;
    struct file_range whole_file;
    int n_models = 0, n_chains = 0, j0, n_new_chains;
    freesasa_structure **ss = NULL, **ssb, *s;

    if( ! (options & FREESASA_SEPARATE_MODELS ||
           options & FREESASA_SEPARATE_CHAINS) ) {
        fail_msg("Options need to specify at least one of FREESASA_SEPARATE_CHAINS "
                 "and FREESASA_SEPARATE_MODELS.");
        return NULL;
    }

    whole_file = freesasa_whole_file(pdb);
    n_models = freesasa_pdb_get_models(pdb,&models);

    if (n_models == FREESASA_FAIL) {
        fail_msg("Problems reading PDB-file.");
        return NULL;
    }
    if (n_models == 0) {
        models = &whole_file;
        n_models = 1;
    }

    //only keep first model if option not provided
    if (! (options & FREESASA_SEPARATE_MODELS) ) n_models = 1;

    //for each model read chains if requested
    if (options & FREESASA_SEPARATE_CHAINS) {
        for (int i = 0; i < n_models; ++i) {
            chains = NULL;
            n_new_chains = freesasa_pdb_get_chains(pdb, models[i], &chains, options);

            if (n_new_chains == FREESASA_FAIL) goto cleanup;
            if (n_new_chains == 0) {
                freesasa_warn("in %s(): No chains found (in model %d).", __func__, i+1);
                continue;
            }

            ssb = ss;
            ss = realloc(ss,sizeof(freesasa_structure*)*(n_chains + n_new_chains));

            if (!ss) {
                ss = ssb;
                mem_fail();
                goto cleanup;
            }

            j0 = n_chains;
            n_chains += n_new_chains;

            for (int j = 0; j < n_new_chains; ++j) ss[j0+j] = NULL;

            for (int j = 0; j < n_new_chains; ++j) {
                s = from_pdb_impl(pdb, chains[j], classifier, options);
                if (s == NULL) goto cleanup;
                ss[j0+j] = s;
                ss[j0+j]->model = ss[n_chains-n_new_chains]->model;
            }

            free(chains);
        }
        *n = n_chains;
    } else {
        ss = malloc(sizeof(freesasa_structure*)*n_models);
        if (!ss) {
            mem_fail();
            goto cleanup;
        }

        for (int i = 0; i < n_models; ++i) ss[i] = NULL;
        *n = n_models;

        for (int i = 0; i < n_models; ++i) {
            s = from_pdb_impl(pdb, models[i], classifier, options);
            if (s == NULL) goto cleanup;
            ss[i] = s;
        }
    }

    if (*n == 0) goto cleanup;

    if (models != &whole_file) free(models);

    return ss;

 cleanup:
    if (ss) for (int i = 0; i < *n; ++i) freesasa_structure_free(ss[i]);
    if (models != &whole_file) free(models);
    free(chains);
    *n = 0;
    free(ss);
    return NULL;
}

freesasa_structure*
freesasa_structure_get_chains(const freesasa_structure *s,
                              const char* chains)
{
    assert(s);
    if (strlen(chains) == 0) return NULL;

    freesasa_structure *new_s = freesasa_structure_new();
    
    if (!new_s) {
        mem_fail();
        return NULL;
    }
    
    new_s->model = s->model;

    for (int i = 0; i < s->number_atoms; ++i) {
        struct atom *ai = s->a[i];
        char c = ai->chain_label;
        if (strchr(chains,c) != NULL) {
            const double *v = freesasa_coord_i(s->xyz,i);
            int res = freesasa_structure_add_atom(new_s, ai->atom_name,
                                                  ai->res_name, ai->res_number,
                                                  c, v[0], v[1], v[2]);
            if (res == FREESASA_FAIL) {
                fail_msg("");
                freesasa_structure_free(new_s);
                return NULL;
            }
        }
    }

    if (new_s->number_atoms == 0) {
        freesasa_structure_free(new_s);
        new_s = NULL;
    }

    return new_s;
}

const char *
freesasa_structure_chain_labels(const freesasa_structure *structure)
{
    assert(structure);
    return structure->chains;
}

const coord_t *
freesasa_structure_xyz(const freesasa_structure *s)
{
    assert(s);
    return s->xyz;
}

int
freesasa_structure_n(const freesasa_structure *s)
{
    assert(s);
    return s->number_atoms;
}

int
freesasa_structure_n_residues(const freesasa_structure *s)
{
    assert(s);
    return s->number_residues;
}

const char *
freesasa_structure_atom_name(const freesasa_structure *s,
                             int i)
{
    assert(s);
    assert(i < s->number_atoms && i >= 0);
    return s->a[i]->atom_name;
}

const char*
freesasa_structure_atom_res_name(const freesasa_structure *s,
                                 int i)
{
    assert(s);
    assert(i < s->number_atoms && i >= 0);
    return s->a[i]->res_name;
}

const char*
freesasa_structure_atom_res_number(const freesasa_structure *s,
                                   int i)
{
    assert(s);
    assert(i < s->number_atoms && i >= 0);
    return s->a[i]->res_number;
}

char
freesasa_structure_atom_chain(const freesasa_structure *s,
                              int i)
{
    assert(s);
    assert(i < s->number_atoms && i >= 0);
    return s->a[i]->chain_label;
}
const char*
freesasa_structure_atom_symbol(const freesasa_structure *s,
                               int i)
{
    assert(s);
    assert(i < s->number_atoms && i >= 0);
    return s->a[i]->symbol;
}

const char*
freesasa_structure_atom_descriptor(const freesasa_structure *s,
                                   int i)
{
    assert(s);
    assert(i < s->number_atoms && i >= 0);
    return s->a[i]->descriptor;
}

int
freesasa_structure_residue_atoms(const freesasa_structure *s,
                                 int r_i,
                                 int *first,
                                 int *last)
{
    assert(s); assert(first); assert(last);
    const int naa = s->number_residues;
    assert(r_i >= 0 && r_i < naa);

    *first = s->res_first_atom[r_i];
    if (r_i == naa-1) *last = s->number_atoms-1;
    else *last = s->res_first_atom[r_i+1]-1;
    assert(*last >= *first);

    return FREESASA_SUCCESS;
}

const char*
freesasa_structure_residue_descriptor(const freesasa_structure *s,
                                      int r_i)
{
    assert(s);
    assert(r_i < s->number_residues && r_i >= 0);
    return s->res_desc[r_i];
}

const char*
freesasa_structure_residue_name(const freesasa_structure *s,
                                int r_i)
{
    assert(s);
    assert(r_i < s->number_residues && r_i >= 0);
    return s->a[s->res_first_atom[r_i]]->res_name;
}

const char*
freesasa_structure_residue_number(const freesasa_structure *s,
                                  int r_i)
{
    assert(s);
    assert(r_i < s->number_residues && r_i >= 0);
    return s->a[s->res_first_atom[r_i]]->res_number;
}

char
freesasa_structure_residue_chain(const freesasa_structure *s,
                                 int r_i)
{
    assert(s);
    assert(r_i < s->number_residues && r_i >= 0);

    return s->a[s->res_first_atom[r_i]]->chain_label;
}

int
freesasa_structure_n_chains(const freesasa_structure *s)
{
    return s->number_chains;
}

int
freesasa_structure_chain_index(const freesasa_structure *s,
                               char chain)
{
    for (int i = 0; i < s->number_chains; ++i) {
        if (s->chains[i] == chain) return i;
    }
    return freesasa_fail("in %s: Chain %c not found.", __func__, chain);
}

int
freesasa_structure_chain_atoms(const freesasa_structure *s,
                               char chain,
                               int *first,
                               int *last)
{
    int c_i = freesasa_structure_chain_index(s,chain),
        n = freesasa_structure_n_chains(s);
    if (c_i < 0) return fail_msg("");

    *first = s->chain_first_atom[c_i];
    if (c_i == n - 1) *last = s->number_atoms-1;
    else *last = s->chain_first_atom[c_i+1] - 1;
    assert(*last >= *first);

    return FREESASA_SUCCESS;
}

int
freesasa_write_pdb(FILE *output,
                   freesasa_result *result,
                   const freesasa_structure *s)
{
    assert(s);
    assert(output);
    assert(result);
    assert(result->sasa);

    const double* values = result->sasa;
    const double* radii = s->radius;
    char buf[PDB_LINE_STRL+1], buf2[6];
    int n = freesasa_structure_n(s);
    if (s->model > 0) fprintf(output,"MODEL     %4d\n",s->model);
    else fprintf(output,             "MODEL        1\n");

    // Write ATOM entries
    for (int i = 0; i < n; ++i) {
        if (s->a[i]->line == NULL) {
            return freesasa_fail("in %s(): PDB input not valid or not present.",
                                 __func__);
        }
        strncpy(buf,s->a[i]->line,PDB_LINE_STRL);
        sprintf(&buf[54],"%6.2f%6.2f",radii[i],values[i]);
        fprintf(output,"%s\n",buf);
    }

    // Write TER  and ENDMDL lines
    errno = 0;
    strncpy(buf2,&buf[6],5);
    buf2[5]='\0';
    fprintf(output,"TER   %5d     %4s %c%4s\nENDMDL\n",
            atoi(buf2)+1, s->a[n-1]->res_name,
            s->a[n-1]->chain_label, s->a[n-1]->res_number);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int
freesasa_structure_model(const freesasa_structure *structure)
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
    return structure->radius;
}

void
freesasa_structure_set_radius(freesasa_structure *structure, const double* radii)
{
    assert(structure);
    assert(radii);
    memcpy(structure->radius, radii, structure->number_atoms*sizeof(double));
}
