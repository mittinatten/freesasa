#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdlib.h>
#if HAVE_STRINGS_H
#include <strings.h>
#endif
#include <errno.h>

#include "classifier.h"
#include "freesasa_internal.h"
#include "pdb.h"

#define STD_CLASSIFIER_NAME "no-name-given"

#define MAX_LINE_LEN 256

/**
    In this file the concept class refers to polar/apolar and type to
    aliphatic/aromatic/etc. See the example configurations in share/.
 */

static const struct classifier_types empty_types = {0, NULL, NULL, NULL};

static const struct classifier_residue empty_residue = {0, NULL, NULL, NULL, NULL, {NULL, 0, 0, 0, 0, 0}};

static const struct freesasa_classifier empty_config = {0, NULL, NULL, NULL};

struct classifier_types *
freesasa_classifier_types_new(void)
{
    struct classifier_types *t = malloc(sizeof(struct classifier_types));
    if (t == NULL)
        mem_fail();
    else
        *t = empty_types;
    return t;
}

void freesasa_classifier_types_free(struct classifier_types *t)
{
    int i;

    if (t != NULL) {
        free(t->type_radius);
        free(t->type_class);
        if (t->name)
            for (i = 0; i < t->n_types; ++i)
                free(t->name[i]);
        free(t->name);

        free(t);
    }
}

struct classifier_residue *
freesasa_classifier_residue_new(const char *name)
{
    struct classifier_residue *res;
    assert(strlen(name) > 0);

    res = malloc(sizeof(struct classifier_residue));

    if (res == NULL)
        mem_fail();
    else {
        *res = empty_residue;
        res->name = strdup(name);
        if (res->name == NULL) {
            mem_fail();
            free(res);
            res = NULL;
        }
    }
    return res;
}

void freesasa_classifier_residue_free(struct classifier_residue *res)
{
    int i;
    if (res != NULL) {
        free(res->name);

        if (res->atom_name)
            for (i = 0; i < res->n_atoms; ++i)
                free(res->atom_name[i]);
        free(res->atom_name);

        free(res->atom_radius);
        free(res->atom_class);

        free(res);
    }
}

freesasa_classifier *
freesasa_classifier_new()
{
    struct freesasa_classifier *cfg = malloc(sizeof(struct freesasa_classifier));
    if (cfg == NULL)
        mem_fail();
    else
        *cfg = empty_config;
    return cfg;
}

void freesasa_classifier_free(freesasa_classifier *c)
{
    int i;
    if (c != NULL) {
        if (c->residue)
            for (i = 0; i < c->n_residues; ++i)
                freesasa_classifier_residue_free(c->residue[i]);
        free(c->residue);
        free(c->residue_name);
        free(c->name);
        free(c);
    }
}

/* check if array of strings has a string that matches key,
   ignores trailing and leading whitespace */
static int
find_string(char **array,
            const char *key,
            int array_size)
{
    int n, i, found = 0;
    char *key_trimmed;

    if (array == NULL || array_size == 0) return -1;

    n = strlen(key);
    key_trimmed = malloc(n + 1);

    if (key_trimmed == NULL) return mem_fail();

    /* remove trailing and leading whitespace */
    sscanf(key, "%s", key_trimmed);

    for (i = 0; i < array_size; ++i) {
        assert(array[i]);
        if (strcmp(array[i], key_trimmed) == 0) {
            found = 1;
            break;
        }
    }

    free(key_trimmed);

    if (found) return i;
    return FREESASA_FAIL;
}

/**
   Removes comments and strips leading and trailing
   whitespace. Returns the length of the stripped line on success,
   FREESASA_FAIL if malloc/realloc fails. Result will be stored in
   the string line, which is assumed to have size MAX_LINE_LEN + 1.
 */
static int
strip_line(char *line,
           const char *input)
{
    char *comment, *first, *last;
    char linebuf[MAX_LINE_LEN + 1];

    assert(strlen(input) <= MAX_LINE_LEN);

    strcpy(linebuf, input);
    comment = strchr(linebuf, '#');
    if (comment) *comment = '\0'; /* skip comments */

    first = linebuf;
    last = linebuf + strlen(linebuf) - 1;
    while (*first == ' ' || *first == '\t')
        ++first;

    if (last > first)
        while (*last == ' ' || *last == '\t' || *last == '\n')
            --last;

    if (first >= last) {
        line[0] = '\0';
        return 0;
    }

    *(last + 1) = '\0';
    strncpy(line, first, MAX_LINE_LEN);

    return strlen(line);
}

/**
    Essentially a safer fscanf(input, "%s", str) limited to the
    current line in input.  Stores the result in 'str' (which should
    be able to store a string of length MAX_LINE_LEN)
*/
static int
get_next_string(FILE *input, char *str)
{
    char line[MAX_LINE_LEN + 1];
    long pos = ftell(input);

    if (fgets(line, MAX_LINE_LEN + 1, input) == NULL) {
        if (ferror(input)) {
            return freesasa_fail(strerror(errno));
        }
        return 0;
    }

    str[0] = '\0';

    sscanf(line, "%s", str);

    fseek(input, pos + strlen(str), SEEK_SET);
    return strlen(str);
}

/**
    Allocates space and stores a line stripped of comments in the line
    pointer. Returns the length of the line on success, FREESASA_FAIL
    for I/O errors.
 */
static int
next_line(char *line,
          FILE *fp)
{
    char linebuf[MAX_LINE_LEN + 1];

    if (fgets(linebuf, MAX_LINE_LEN + 1, fp) == NULL) {
        if (ferror(fp)) {
            return fail_msg(strerror(errno));
        }

        if (feof(fp)) {
            line[0] = '\0';
            return 0;
        }
    }

    return strip_line(line, linebuf);
}

/**
    Find offset of str in line, returns -1 if not found
    Ignores comments
*/
static inline int
locate_string(const char *line,
              const char *str)
{
    int NOT_FOUND = -1;
    char *loc, buf[MAX_LINE_LEN + 1];

    assert(line);
    assert(strlen(line) <= MAX_LINE_LEN);
    assert(str);

    if (strlen(line) == 0) {
        return NOT_FOUND;
    }

    strcpy(buf, line);

    /* skip comments */
    loc = strstr(buf, "#");
    if (loc == buf) {
        return NOT_FOUND;
    } else if (loc != NULL) {
        *loc = '\0';
    }

    loc = strstr(buf, str);
    if (loc != NULL) {
        return loc - buf;
    }

    return NOT_FOUND;
}

/**
    If string exists on line its location is stored in this_range, and
    if prev_range is non-null it is set to end at the same location.
 */
static inline int
try_register_stringloc(const char *line,
                       const char *str,
                       long last_tell,
                       struct file_range *this_range,
                       struct file_range **prev_range)
{
    int pos, NOT_FOUND = -1;

    if (strlen(line) == 0) return NOT_FOUND;

    pos = locate_string(line, str);
    if (pos >= 0) {
        this_range->begin = last_tell + pos;
        if (*prev_range) (*prev_range)->end = last_tell + pos;
        (*prev_range) = this_range;
        return last_tell + pos;
    }
    return NOT_FOUND;
}

/**
    Checks that input file has the required fields and locates the
    'types' and 'atoms' sections. No syntax checking. Return
    FREESASA_SUCCESS if file seems ok, FREESASA_FAIL if either/both of
    the sections are missing, or file invalid in some other way, or there was an error reading the file.
 */
static int
check_file(FILE *input,
           struct file_range *types,
           struct file_range *atoms,
           struct file_range *name)
{
    long last_tell;
    char line[MAX_LINE_LEN + 1];
    struct file_range *last_range = NULL;

    assert(input);
    assert(types);
    assert(atoms);

    last_tell = ftell(input);

    /* this allows us to detect wether a section wasn't found later */
    types->begin = atoms->begin = name->begin = -1;

    while (fgets(line, MAX_LINE_LEN + 1, input)) {
        try_register_stringloc(line, "types:", last_tell, types, &last_range);
        try_register_stringloc(line, "atoms:", last_tell, atoms, &last_range);
        try_register_stringloc(line, "name:", last_tell, name, &last_range);
        last_tell = ftell(input);

        if (strlen(line) == MAX_LINE_LEN &&
            line[MAX_LINE_LEN - 1] != '\n') {
            return fail_msg("Lines in classifier files can only be %d characters or less",
                            MAX_LINE_LEN);
        }
    }

    if (ferror(input)) {
        return fail_msg(strerror(errno));
    }

    if (last_range != NULL) {
        last_range->end = last_tell;
    }
    rewind(input);

    if (name->begin == -1) {
        freesasa_warn("input configuration lacks the entry 'name:', "
                      "will use '" STD_CLASSIFIER_NAME "'");
    }

    if ((types->begin == -1) ||
        (atoms->begin == -1)) {
        return fail_msg("input configuration lacks (at least) one of "
                        "the entries 'types:' or 'atoms:'");
    }

    return FREESASA_SUCCESS;
}

int freesasa_classifier_parse_class(const char *name)
{
#if HAVE_STRNCASECMP
    if (strncasecmp(name, "apolar", 6) == 0) {
        return FREESASA_ATOM_APOLAR;
    } else if (strncasecmp(name, "polar", 5) == 0) {
        return FREESASA_ATOM_POLAR;
    } else {
        return fail_msg("only atom classes allowed are 'polar' and 'apolar'"
                        " (case insensitive)");
    }
#else
    if (strncmp(name, "apolar", 6) == 0) {
        return FREESASA_ATOM_APOLAR;
    } else if (strncmp(name, "polar", 5) == 0) {
        return FREESASA_ATOM_POLAR;
    } else {
        return fail_msg("only atom classes allowed are 'polar' and 'apolar'");
    }
#endif
}

/**
    Add type. Returns the index of the new type on success,
    FREESASA_FAIL if realloc/strdup fails, FREESASA_WARN if type
    already known (ignore duplicates).
 */
int freesasa_classifier_add_type(struct classifier_types *types,
                                 const char *type_name,
                                 const char *class_name,
                                 double r)
{
    int the_class;
    int n = types->n_types + 1;
    char **tn = types->name;
    double *tr = types->type_radius;
    freesasa_atom_class *tc = types->type_class;

    if (find_string(types->name, type_name, types->n_types) >= 0)
        return freesasa_warn("ignoring duplicate configuration entry for '%s'", type_name);

    the_class = freesasa_classifier_parse_class(class_name);
    if (the_class == FREESASA_FAIL) return fail_msg("");

    if ((types->name = realloc(tn, sizeof(char *) * n)) == NULL) {
        types->name = tn;
        return mem_fail();
    }

    if ((types->type_radius = realloc(tr, sizeof(double) * n)) == NULL) {
        types->type_radius = tr;
        return mem_fail();
    }

    if ((types->type_class = realloc(tc, sizeof(int) * n)) == NULL) {
        types->type_class = tc;
        return mem_fail();
    }

    if ((types->name[n - 1] = strdup(type_name)) == NULL) {
        return mem_fail();
    }

    types->n_types++;
    types->type_radius[types->n_types - 1] = r;
    types->type_class[types->n_types - 1] = the_class;

    return types->n_types - 1;
}

/**
    Read a line specifying a type, store it in the config. Returns
    warning for duplicates, failures for syntax errors or memory
    allocation errors.
 */
static int
read_types_line(struct classifier_types *types,
                const char *line)
{
    int the_type, ret = FREESASA_SUCCESS;
    double r;
    char buf1[MAX_LINE_LEN + 1], buf2[MAX_LINE_LEN + 1];

    assert(strlen(line) <= MAX_LINE_LEN);

    if (sscanf(line, "%s %lf %s", buf1, &r, buf2) == 3) {
        the_type = freesasa_classifier_add_type(types, buf1, buf2, r);
        if (the_type == FREESASA_FAIL) ret = fail_msg("");
        if (the_type == FREESASA_WARN) ret = FREESASA_WARN;
    } else {
        ret = fail_msg("could not parse line '%s' in configuration, "
                       "expecting triplet of type 'TYPE [RADIUS] CLASS' for "
                       "example 'C_ALI 2.00 apolar'",
                       line);
    }

    return ret;
}

/**
    Reads info about types from the user config. Associates each type
    with a class and a radius in the config struct. Returns
    FREESASA_SUCCESS on success, FREESASA_FAIL on syntax or memory
    allocation errors.
 */
static int
read_types(struct classifier_types *types,
           FILE *input,
           struct file_range fi)
{
    char line[MAX_LINE_LEN + 1];
    int ret = FREESASA_SUCCESS, nl;

    fseek(input, fi.begin, SEEK_SET);

    /* read command (and discard) */
    if (next_line(line, input) > 0) {
        char buf[7]; /* we should not get here if the line isn't "types:" (plus whitespace) */
        if (sscanf(line, "%6s", buf) == 0) return FREESASA_FAIL;
        assert(strcmp(buf, "types:") == 0);
    } else {
        return FREESASA_FAIL;
    }

    while (ftell(input) < fi.end) {
        nl = next_line(line, input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) {
            ret = nl;
            break;
        };
        ret = read_types_line(types, line);
        if (ret == FREESASA_FAIL) break;
    }

    return ret;
}

/**
    Add atom to residue. Returns index of the new atom on
    success. FREESASA_FAIL if memory allocation fails. FREESASA_WARN
    if the atom has already been added.
 */
int freesasa_classifier_add_atom(struct classifier_residue *res,
                                 const char *name,
                                 double radius,
                                 int the_class)
{
    int n;
    char **an = res->atom_name;
    double *ar = res->atom_radius;
    freesasa_atom_class *ac = res->atom_class;

    if (find_string(res->atom_name, name, res->n_atoms) >= 0)
        return freesasa_warn("ignoring duplicate configuration entry for atom '%s %s'",
                             res->name, name);
    n = res->n_atoms + 1;

    if ((res->atom_name = realloc(res->atom_name, sizeof(char *) * n)) == NULL) {
        res->atom_name = an;
        return mem_fail();
    }
    if ((res->atom_radius = realloc(res->atom_radius, sizeof(double) * n)) == NULL) {
        res->atom_radius = ar;
        return mem_fail();
    }
    if ((res->atom_class = realloc(res->atom_class, sizeof(int) * n)) == NULL) {
        res->atom_class = ac;
        return mem_fail();
    }
    if ((res->atom_name[n - 1] = strdup(name)) == NULL)
        return mem_fail();

    ++res->n_atoms;
    res->atom_radius[n - 1] = radius;
    res->atom_class[n - 1] = the_class;

    return n - 1;
}

/**
    Add residue to config. If the residue already exists, it returns
    the index of that residue, else it returns the index of the new
    residue. Returns FREESASA_FAILURE if realloc/strdup fails.
 */
int freesasa_classifier_add_residue(struct freesasa_classifier *c,
                                    const char *name)
{
    char **rn = c->residue_name;
    struct classifier_residue **cr = c->residue;
    int res = find_string(c->residue_name, name, c->n_residues);

    if (res >= 0) return res;

    res = c->n_residues + 1;

    if ((c->residue_name = realloc(rn, sizeof(char *) * res)) == NULL) {
        c->residue_name = rn;
        return mem_fail();
    }

    if ((c->residue = realloc(cr, sizeof(struct classifier_residue *) * res)) == NULL) {
        c->residue = cr;
        return mem_fail();
    }

    if ((c->residue[res - 1] = freesasa_classifier_residue_new(name)) == NULL) {
        return mem_fail();
    }

    ++c->n_residues;
    c->residue_name[res - 1] = c->residue[res - 1]->name;
    return res - 1;
}

/**
    Read a line specifying an atom, store it in the config. Use
    supplied types to add assign radius and class. Returns
    FREESASA_WARN for duplicates. Returns FREESASA_FAIL for syntax
    errors or memory allocation errors. FREESASA_SUCCESS else.
 */
static int
read_atoms_line(struct freesasa_classifier *c,
                const struct classifier_types *types,
                const char *line)
{
    char buf1[MAX_LINE_LEN + 1], buf2[MAX_LINE_LEN + 1], buf3[MAX_LINE_LEN + 1];
    int res, type, atom;

    assert(strlen(line) <= MAX_LINE_LEN);

    if (sscanf(line, "%s %s %s", buf1, buf2, buf3) == 3) {
        if (strlen(buf1) > PDB_ATOM_RES_NAME_STRL) {
            return fail_msg("residue name %s is too long in classifier file", buf1);
        }

        if (strlen(buf2) > PDB_ATOM_NAME_STRL) {
            return fail_msg("atom name %s is too long in classifier file", buf2);
        }

        type = find_string(types->name, buf3, types->n_types);

        if (type < 0) {
            return fail_msg("unknown atom type '%s' in configuration, line '%s'",
                            buf3, line);
        }

        res = freesasa_classifier_add_residue(c, buf1);

        if (res == FREESASA_FAIL) return fail_msg("");

        atom = freesasa_classifier_add_atom(c->residue[res],
                                            buf2,
                                            types->type_radius[type],
                                            types->type_class[type]);

        if (atom == FREESASA_FAIL) return fail_msg("");
        if (atom == FREESASA_WARN) return FREESASA_WARN;

    } else {
        return fail_msg("could not parse configuration, line '%s', "
                        "expecting triplet of type "
                        "'RESIDUE ATOM CLASS', for example 'ALA CB C_ALI'",
                        line);
    }
    return FREESASA_SUCCESS;
}

/**
    Reads atom configurations from config-file. Associates each atom
    with a radius and class using the types that should already have
    been stored in the config struct.
 */
static int
read_atoms(struct freesasa_classifier *c,
           struct classifier_types *types,
           FILE *input,
           struct file_range fi)
{
    char line[MAX_LINE_LEN + 1], buf[MAX_LINE_LEN + 1];
    int ret = FREESASA_SUCCESS, nl;

    fseek(input, fi.begin, SEEK_SET);

    /* read command (and discard) */
    if (next_line(line, input) > 0) {
        assert(strlen(line) <= MAX_LINE_LEN);
        if (sscanf(line, "%s", buf) == 0) return FREESASA_FAIL;
        assert(strcmp(buf, "atoms:") == 0);
    } else {
        return FREESASA_FAIL;
    }

    while (ftell(input) < fi.end) {
        nl = next_line(line, input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) return fail_msg("");
        ret = read_atoms_line(c, types, line);
        if (ret == FREESASA_FAIL) break;
    }

    return ret;
}

static int
read_name(struct freesasa_classifier *classifier,
          FILE *input,
          struct file_range fi)
{
    char buf[MAX_LINE_LEN + 1];

    if (fi.begin < 0)
        return FREESASA_SUCCESS; /* name not set? */

    fseek(input, fi.begin, SEEK_SET);
    if (get_next_string(input, buf) <= 0)
        return fail_msg("");

    assert(strcmp(buf, "name:") == 0);

    if (get_next_string(input, buf) <= 0) {
        return fail_msg("empty name for configuration?");
    }

    classifier->name = strdup(buf);
    if (classifier->name == NULL) {
        return mem_fail();
    }

    return FREESASA_SUCCESS;
}

static struct freesasa_classifier *
read_config(FILE *input)
{
    struct file_range types_section, atoms_section, name_section;
    struct freesasa_classifier *classifier = NULL;
    struct classifier_types *types = NULL;

    assert(input);

    if (!(types = freesasa_classifier_types_new()))
        goto cleanup;
    if (!(classifier = freesasa_classifier_new()))
        goto cleanup;
    if (check_file(input, &types_section, &atoms_section, &name_section))
        goto cleanup;
    if (read_name(classifier, input, name_section))
        goto cleanup;
    if (read_types(types, input, types_section))
        goto cleanup;
    if (read_atoms(classifier, types, input, atoms_section))
        goto cleanup;

    freesasa_classifier_types_free(types);

    return classifier;

cleanup:
    freesasa_classifier_free(classifier);
    freesasa_classifier_types_free(types);
    return NULL;
}

/**
    See if an atom_name has been defined for the residue ANY (writes
    indices to the provided pointers).
 */
static void
find_any(const struct freesasa_classifier *c,
         const char *atom_name,
         int *res, int *atom)
{
    *res = find_string(c->residue_name, "ANY", c->n_residues);
    if (*res >= 0) {
        *atom = find_string(c->residue[*res]->atom_name,
                            atom_name,
                            c->residue[*res]->n_atoms);
    }
}
/**
    Find the residue and atom index of an atom in the supplied
    configuration. Prints error and returns FREESASA_WARN if not
    found.
 */
static int
find_atom(const struct freesasa_classifier *c,
          const char *res_name,
          const char *atom_name,
          int *res,
          int *atom)
{
    const struct classifier_residue *residue;
    *atom = -1;
    *res = find_string(c->residue_name, res_name, c->n_residues);
    if (*res < 0) {
        find_any(c, atom_name, res, atom);
    } else {
        residue = c->residue[*res];
        *atom = find_string(residue->atom_name, atom_name, residue->n_atoms);
        if (*atom < 0) {
            find_any(c, atom_name, res, atom);
        }
    }
    if (*atom < 0) {
        return FREESASA_WARN;
    }
    return FREESASA_SUCCESS;
}

double
freesasa_classifier_radius(const freesasa_classifier *classifier,
                           const char *res_name,
                           const char *atom_name)
{
    int res, atom, status;

    assert(classifier);
    assert(res_name);
    assert(atom_name);

    status = find_atom(classifier, res_name, atom_name, &res, &atom);
    if (status == FREESASA_SUCCESS)
        return classifier->residue[res]->atom_radius[atom];
    return -1.0;
}

freesasa_atom_class
freesasa_classifier_class(const freesasa_classifier *classifier,
                          const char *res_name,
                          const char *atom_name)
{
    int res, atom, status;

    assert(classifier);
    assert(res_name);
    assert(atom_name);

    status = find_atom(classifier, res_name, atom_name, &res, &atom);
    if (status == FREESASA_SUCCESS)
        return classifier->residue[res]->atom_class[atom];
    return FREESASA_ATOM_UNKNOWN;
}

const char *
freesasa_classifier_class2str(freesasa_atom_class atom_class)
{
    switch (atom_class) {
    case FREESASA_ATOM_APOLAR:
        return "Apolar";
    case FREESASA_ATOM_POLAR:
        return "Polar";
    case FREESASA_ATOM_UNKNOWN:
        return "Unknown";
    }
    fail_msg("invalid atom class");
    return NULL;
}

freesasa_nodearea
freesasa_result_classes(const freesasa_structure *structure,
                        const freesasa_result *result)
{
    freesasa_nodearea area = {"whole-structure", 0, 0, 0, 0, 0};
    freesasa_range_nodearea(&area, structure, result,
                            0, freesasa_structure_n(structure) - 1);
    return area;
}

freesasa_classifier *
freesasa_classifier_from_file(FILE *file)
{
    struct freesasa_classifier *classifier = read_config(file);

    if (classifier == NULL) {
        fail_msg("");
        return NULL;
    }

    return classifier;
}

const freesasa_nodearea *
freesasa_classifier_residue_reference(const freesasa_classifier *classifier,
                                      const char *res_name)
{
    int res = find_string(classifier->residue_name, res_name, classifier->n_residues);

    if (res < 0) return NULL;

    return &classifier->residue[res]->max_area;
}

const char *
freesasa_classifier_name(const freesasa_classifier *classifier)
{
    return classifier->name;
}

struct symbol_radius {
    const char symbol[3];
    double radius;
};

/* Taken from:

   Mantina et al. "Consistent van der Waals Radii for
   the Whole Main Group". J. Phys. Chem. A, 2009, 113 (19), pp
   5806–5812.

   Many of these elements, if they occur in a PDB file, should
   probably rather be skipped than used in a SASA calculation, and
   ionization will change the effective radius.

*/
static const struct symbol_radius symbol_radius[] = {
    /* elements that actually occur in the regular amino acids and nucleotides */
    {" H", 1.10},
    {" C", 1.70},
    {" N", 1.55},
    {" O", 1.52},
    {" P", 1.80},
    {" S", 1.80},
    {"SE", 1.90},
    /* some others, values pulled from gemmi elem.hpp */
    /* Halogens */
    {" F", 1.47},
    {"CL", 1.75},
    {"BR", 1.83},
    {" I", 1.98},
    /* Alkali and Alkali Earth metals */
    {"LI", 1.81},
    {"BE", 1.53},
    {"NA", 2.27},
    {"MG", 1.73},
    {" K", 2.75},
    {"CA", 2.31},
    {"RB", 3.03},
    {"SR", 2.49},
    {"CS", 3.43},
    {"BA", 2.68},
    {"FR", 3.48},
    {"RA", 2.83},
    /* Transition metals */
    {"SC", 2.11},
    {"TI", 1.95},
    {" V", 1.06},
    {"CR", 1.13},
    {"MN", 1.19},
    {"FE", 1.26},
    {"CO", 1.13},
    {"NI", 1.63},
    {"CU", 1.40},
    {"ZN", 1.39},
    {" Y", 1.61},
    {"ZR", 1.42},
    {"NB", 1.33},
    {"MO", 1.75},
    {"TC", 2.00},
    {"RU", 1.20},
    {"RH", 1.22},
    {"PD", 1.63},
    {"AG", 1.72},
    {"CD", 1.58},
    {"HF", 1.40},
    {"TA", 1.22},
    {" W", 1.26},
    {"RE", 1.30},
    {"OS", 1.58},
    {"IR", 1.22},
    {"PT", 1.75},
    {"AU", 1.66},
    {"HG", 1.55},
    /* Post-Transition metals */
    {"AL", 1.84},
    {"GA", 1.87},
    {"IN", 1.93},
    {"SN", 2.17},
    {"TL", 1.96},
    {"PB", 2.02},
    {"BI", 2.07},
    {"PO", 1.97},
    /* Metalloid */
    {" B", 1.92},
    {"SI", 2.10},
    {"GE", 2.11},
    {"AS", 1.85},
    {"SB", 2.06},
    {"TE", 2.06},
    {"AT", 2.02},
    /* Noble gases */
    {"HE", 1.40},
    {"NE", 1.54},
    {"AR", 1.88},
    {"KR", 2.02},
    {"XE", 2.16},
    {"RN", 2.20},
    /* Lanthanoids */
    {"LA", 1.83},
    {"CE", 1.86},
    {"PR", 1.62},
    {"ND", 1.79},
    {"PM", 1.76},
    {"SM", 1.74},
    {"EU", 1.96},
    {"GD", 1.69},
    {"TB", 1.66},
    {"DY", 1.63},
    {"HO", 1.61},
    {"ER", 1.59},
    {"TM", 1.57},
    {"YB", 1.54},
    {"LU", 1.53},
    /* Actinoids */
    {"AC", 2.12},
    {"TH", 1.84},
    {"PA", 1.60},
    {" U", 1.86},
    {"NP", 1.71},
    {"PU", 1.67},
    {"AM", 1.66},
    {"CM", 1.65},
    {"BK", 1.64},
    {"CF", 1.63},
    {"ES", 1.62},
    {"FM", 1.61},
    {"MD", 1.60},
    {"NO", 1.59},
    {"LR", 1.58},
};

double
freesasa_guess_radius(const char *input_symbol)
{
    int n_symbol, i;
    char symbol[3];

    assert(input_symbol);
    snprintf(symbol, 3, "%2s", input_symbol);

    n_symbol = sizeof(symbol_radius) / sizeof(struct symbol_radius);

    for (i = 0; i < n_symbol; ++i) {
        if (strcmp(symbol, symbol_radius[i].symbol) == 0)
            return symbol_radius[i].radius;
    }
    return -1.0;
}

// clang-format off
/* The residue types that are returned by freesasa_classify_residue() */
enum residue {
    /* Regular amino acids */
    ALA = 0, ARG, ASN, ASP,
    CYS, GLN, GLU, GLY,
    HIS, ILE, LEU, LYS,
    MET, PHE, PRO, SER,
    THR, TRP, TYR, VAL,
    /* some non-standard ones */
    CSE, SEC, PYL, PYH, ASX, GLX,
    /* residue unknown */
    RES_UNK,
    /* capping N- and C-terminal groups (usually HETATM) */
    ACE, NH2,
    /* DNA */
    DA, DC, DG, DT, DU, DI,
    /* RNA (avoid one-letter enums) */
    RA, RC, RG, RU, RI, RT,
    /* generic nucleotide */
    NN
};

/* Residue types, make sure this always matches the corresponding enum. */
static const char *residue_names[] = {
    /* amino acids */
    "ALA", "ARG", "ASN", "ASP",
    "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS",
    "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL",
    /* non-standard amino acids */
    "CSE", "SEC", "PYL", "PYH", /* SEC and PYL are standard names, CSE and PYH are found in some early files */
    "ASX", "GLX",
    "UNK",
    /* capping groups */
    "ACE", "NH2",
    /* DNA */
    "DA", "DC", "DG", "DT", "DU", "DI",
    /* RNA */
    "A", "C", "G", "U", "I", "T",
    /* General nucleotide */
    "N"};
// clang-format on

int freesasa_classify_n_residue_types()
{
    return NN + 1;
}

int freesasa_classify_residue(const char *res_name)
{
    int i;
    char cpy[PDB_ATOM_RES_NAME_STRL + 1];

    sscanf(res_name, "%s", cpy);

    for (i = ALA; i < freesasa_classify_n_residue_types(); ++i) {
        if (strcmp(cpy, residue_names[i]) == 0) return i;
    }

    return RES_UNK;
}

const char *
freesasa_classify_residue_name(int residue_type)
{
    assert(residue_type >= 0 && residue_type <= NN);
    return residue_names[residue_type];
}

int freesasa_atom_is_backbone(const char *atom_name)
{
    const char *bb[] = {"CA", "N", "O", "C", "OXT",
                        "P", "OP1", "OP2", "O5'", "C5'", "C4'",
                        "O4'", "C3'", "O3'", "C2'", "C1'"};
    char name[PDB_ATOM_NAME_STRL + 1];
    int i;

    name[0] = '\0';
    sscanf(atom_name, "%s", name); /* trim whitespace */

    if (strlen(name) == 0) return 0;
    for (i = 0; i < sizeof(bb) / sizeof(const char *); ++i) {
        if (strcmp(name, bb[i]) == 0) {
            return 1;
        }
    }
    return 0;
}

#if USE_CHECK
#include <check.h>
#include <math.h>

START_TEST(test_classifier)
{
    struct classifier_types *types = freesasa_classifier_types_new();
    struct classifier_residue *residue_cfg = freesasa_classifier_residue_new("ALA");
    struct freesasa_classifier *clf = freesasa_classifier_new();

    freesasa_set_verbosity(FREESASA_V_SILENT);

    ck_assert_int_eq(freesasa_classifier_parse_class("A"), FREESASA_FAIL);
#if HAVE_STRNCASECMP
    ck_assert_int_eq(freesasa_classifier_parse_class("POLAR"), FREESASA_ATOM_POLAR);
    ck_assert_int_eq(freesasa_classifier_parse_class("APOLAR"), FREESASA_ATOM_APOLAR);
#endif
    ck_assert_int_eq(freesasa_classifier_parse_class("polar"), FREESASA_ATOM_POLAR);
    ck_assert_int_eq(freesasa_classifier_parse_class("apolar"), FREESASA_ATOM_APOLAR);

    ck_assert_int_eq(types->n_types, 0);
    ck_assert_int_eq(freesasa_classifier_add_type(types, "a", "A", 1.0), FREESASA_FAIL);
    ck_assert_int_eq(freesasa_classifier_add_type(types, "a", "polar", 1.0), 0);
    ck_assert_int_eq(freesasa_classifier_add_type(types, "b", "apolar", 2.0), 1);
    ck_assert_int_eq(freesasa_classifier_add_type(types, "b", "polar", 1.0), FREESASA_WARN);
    ck_assert_int_eq(freesasa_classifier_add_type(types, "c", "apolar", 3.0), 2);
    ck_assert_int_eq(types->n_types, 3);
    ck_assert_str_eq(types->name[0], "a");
    ck_assert_str_eq(types->name[1], "b");
    ck_assert_str_eq(types->name[2], "c");
    ck_assert(fabs(types->type_radius[0] - 1.0) < 1e-10);
    ck_assert(fabs(types->type_radius[1] - 2.0) < 1e-10);
    ck_assert(fabs(types->type_radius[2] - 3.0) < 1e-10);

    freesasa_classifier_types_free(types);
    types = freesasa_classifier_types_new();

    ck_assert_int_eq(read_types_line(types, ""), FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types, "a"), FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types, "a 1.0"), FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types, "a b C"), FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types, "a 1.0 C"), FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types, "a 1.0 apolar"), FREESASA_SUCCESS);
    ck_assert_int_eq(read_types_line(types, "b 2.0 polar"), FREESASA_SUCCESS);
    ck_assert_int_eq(types->n_types, 2);
    ck_assert_str_eq(types->name[0], "a");
    ck_assert_str_eq(types->name[1], "b");
    ck_assert(fabs(types->type_radius[0] - 1.0) < 1e-10);
    ck_assert(fabs(types->type_radius[1] - 2.0) < 1e-10);

    ck_assert_int_eq(freesasa_classifier_add_atom(residue_cfg, "C", 1.0, 0), 0);
    ck_assert_int_eq(freesasa_classifier_add_atom(residue_cfg, "CB", 2.0, 0), 1);
    ck_assert_int_eq(freesasa_classifier_add_atom(residue_cfg, "CB", 2.0, 0), FREESASA_WARN);
    ck_assert_str_eq(residue_cfg->atom_name[0], "C");
    ck_assert_str_eq(residue_cfg->atom_name[1], "CB");
    ck_assert(fabs(residue_cfg->atom_radius[0] - 1.0) < 1e-10);
    ck_assert(fabs(residue_cfg->atom_radius[1] - 2.0) < 1e-10);
    freesasa_classifier_residue_free(residue_cfg);

    ck_assert_int_eq(freesasa_classifier_add_residue(clf, "A"), 0);
    ck_assert_int_eq(freesasa_classifier_add_residue(clf, "B"), 1);
    ck_assert_int_eq(freesasa_classifier_add_residue(clf, "B"), 1);
    ck_assert_int_eq(clf->n_residues, 2);
    ck_assert_str_eq(clf->residue_name[0], "A");
    ck_assert_str_eq(clf->residue_name[1], "B");
    ck_assert_str_eq(clf->residue[0]->name, "A");

    freesasa_classifier_free(clf);
    clf = freesasa_classifier_new();

    ck_assert_int_eq(read_atoms_line(clf, types, "A A"), FREESASA_FAIL);
    ck_assert_int_eq(read_atoms_line(clf, types, "A A bla"), FREESASA_FAIL);
    ck_assert_int_eq(read_atoms_line(clf, types, "ALA CA a"), FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf, types, "ALA CB b"), FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf, types, "ARG CA a"), FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf, types, "ARG CB b"), FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf, types, "ARG CG b"), FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf, types, "TOOLONGRESNAME CG b"), FREESASA_FAIL);
    ck_assert_int_eq(read_atoms_line(clf, types, "ARG TOOLONGATOMNAME b"), FREESASA_FAIL);
    ck_assert_int_eq(clf->n_residues, 2);
    ck_assert_str_eq(clf->residue_name[0], "ALA");
    ck_assert_str_eq(clf->residue_name[1], "ARG");
    ck_assert_int_eq(clf->residue[0]->n_atoms, 2);
    ck_assert_str_eq(clf->residue[0]->atom_name[0], "CA");
    ck_assert_str_eq(clf->residue[0]->atom_name[1], "CB");
    ck_assert(fabs(clf->residue[0]->atom_radius[0] - 1.0) < 1e-5);
    ck_assert(fabs(clf->residue[0]->atom_radius[1] - 2.0) < 1e-5);

    freesasa_classifier_free(clf);
    freesasa_classifier_types_free(types);

    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST(test_classifier_utils)
{
    const char *strarr[] = {"A", "B", "C"};
    const char *line[] = {"# Bla", " # Bla", "Bla # Bla", " Bla # Bla", "#Bla #Alb"};

    char dummy_str[MAX_LINE_LEN + 1];
    ck_assert_int_eq(find_string((char **)strarr, "A", 3), 0);
    ck_assert_int_eq(find_string((char **)strarr, "B", 3), 1);
    ck_assert_int_eq(find_string((char **)strarr, "C", 3), 2);
    ck_assert_int_eq(find_string((char **)strarr, "D", 3), -1);
    ck_assert_int_eq(find_string((char **)strarr, " C ", 3), 2);
    ck_assert_int_eq(find_string((char **)strarr, "CC", 3), -1);

    ck_assert_int_eq(strip_line(dummy_str, line[0]), 0);
    ck_assert_int_eq(strip_line(dummy_str, line[1]), 0);
    ck_assert_int_eq(strip_line(dummy_str, line[2]), 3);
    ck_assert_str_eq(dummy_str, "Bla");

    ck_assert_int_eq(strip_line(dummy_str, line[3]), 3);
    ck_assert_str_eq(dummy_str, "Bla");

    ck_assert_int_eq(strip_line(dummy_str, line[4]), 0);

    const char *str = "foo bar # baz";
    ck_assert_int_eq(locate_string(str, "Foo"), -1);
    ck_assert_int_eq(locate_string(str, "foo"), 0);
    ck_assert_int_eq(locate_string(str, "bar"), 4);
    ck_assert_int_eq(locate_string(str, "baz"), -1);
    ck_assert_int_eq(locate_string(str, "ar"), 5);

    struct file_range this_range, *last_range = NULL;
    ck_assert_int_eq(try_register_stringloc(str, "foo", 0, &this_range, &last_range), 0);
    ck_assert_int_eq(this_range.begin, 0);
    ck_assert_ptr_eq(last_range, &this_range);
    ck_assert_int_eq(try_register_stringloc(str, "bar", 0, &this_range, &last_range), 4);
    ck_assert_ptr_eq(last_range, &this_range);
    ck_assert_int_eq(last_range->end, 4);
    ck_assert_int_eq(try_register_stringloc(str, "baz", 0, &this_range, &last_range), -1);
}
END_TEST

TCase *
test_classifier_static()
{
    TCase *tc = tcase_create("classifier.c static");
    tcase_add_test(tc, test_classifier);
    tcase_add_test(tc, test_classifier_utils);

    return tc;
}

#endif /** USE_CHECK */
