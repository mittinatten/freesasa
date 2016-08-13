#if HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#if HAVE_STRINGS_H
#include <strings.h>
#endif
#include <errno.h>
#include "classifier.h"
#include "freesasa_internal.h"

/**
    In this file the concept class refers to polar/apolar and type to
    aliphatic/aromatic/etc. See the example configurations in share/.
 */

static const struct classifier_types empty_types = {0, NULL, NULL, NULL};

static const struct classifier_residue empty_residue = {0, NULL, NULL, NULL, NULL, {NULL, 0, 0, 0, 0, 0}};

static const struct freesasa_classifier empty_config = {0, NULL, NULL, NULL};

struct classifier_types*
freesasa_classifier_types_new()
{
    struct classifier_types *t = malloc(sizeof(struct classifier_types));
    if (t == NULL) mem_fail();
    else *t = empty_types;
    return t;
}

void
freesasa_classifier_types_free(struct classifier_types* t)
{
    if (t != NULL) {
        free(t->type_radius);
        free(t->type_class);
        if (t->name)
            for (int i = 0; i < t->n_types; ++i)
                free(t->name[i]);
        free(t->name);

        free(t);
    }
}

struct classifier_residue*
freesasa_classifier_residue_new(const char* name)
{
    assert(strlen(name) > 0);
    struct classifier_residue *res = malloc(sizeof(struct classifier_residue));
    if (res == NULL) mem_fail();
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

void
freesasa_classifier_residue_free(struct classifier_residue* res)
{
    if (res != NULL) {
        free(res->name);

        if (res->atom_name)
            for (int i = 0; i < res->n_atoms; ++i)
                free(res->atom_name[i]);
        free(res->atom_name);

        free(res->atom_radius);
        free(res->atom_class);

        free(res);
    }
}

freesasa_classifier* 
freesasa_classifier_new()
{
    struct freesasa_classifier *cfg = malloc(sizeof(struct freesasa_classifier));
    if (cfg == NULL) mem_fail();
    else *cfg = empty_config;
    return cfg;
}

void
freesasa_classifier_free(freesasa_classifier *c)
{
    if (c != NULL) {
        if (c->residue)
            for (int i = 0; i < c->n_residues; ++i)
                freesasa_classifier_residue_free(c->residue[i]);
        free(c->residue);
        free(c->residue_name);

        free(c);
    }
}

//! check if array of strings has a string that matches key, ignores trailing and leading whitespace
static int 
find_string(char **array,
            const char *key,
            int array_size)
{
    assert(key);
    if (array == NULL || array_size == 0) return -1;

    int n = strlen(key);
    char key_trimmed[n+1];

    // remove trailing and leading whitespace
    sscanf(key,"%s",key_trimmed);
    for (int i = 0; i < array_size; ++i) {
        assert(array[i]);
        if (strcmp(array[i],key_trimmed) == 0) return i;
    }
    return FREESASA_FAIL;
}

/**
    Checks that input file has the required fields and locates the
    'types' and 'atoms' sections. No syntax checking. Return
    FREESASA_SUCCESS if file seems ok, FREESASA_FILE if either/both of
    the sections are missing.
 */
static int 
check_file(FILE *input,
           struct file_range *types, 
           struct file_range *atoms)
{
    assert(input); assert(types); assert(atoms);
    long last_tell;
    char buf[200];
    struct file_range *last_range = NULL;
        
    last_tell = ftell(input);
    types->begin = atoms->begin = -1;
    while (fscanf(input,"%s",buf) > 0) {
        if (strcmp(buf,"types:") == 0) {
            types->begin = last_tell;
            if (last_range) last_range->end = last_tell;
            last_range = types;
        }
        if (strcmp(buf,"atoms:") == 0) {
            atoms->begin = last_tell;
            if (last_range) last_range->end = last_tell;
            last_range = atoms;
        }
        last_tell = ftell(input);
    }
    if (last_range != NULL) { 
        last_range->end = last_tell;
    } 
    rewind(input);

    if ((types->begin == -1) || 
        (atoms->begin == -1)) {
        return fail_msg("Input configuration lacks (at least) one of "
                        "the entries 'types:' or 'atoms:'.");
    }
    
    return FREESASA_SUCCESS;
}

/**
   Removes comments and strips leading and trailing
   whitespace. Returns the length of the stripped line on success,
   FREESASA_FAIL if malloc/realloc fails.
 */
static int
strip_line(char **line,
           const char *input) 
{
    int n = strlen(input);
    char linebuf[n+1];
    char *comment, *first, *last, *tmp;
    
    strncpy(linebuf, input, strlen(input)+1);
    comment = strchr(linebuf, '#');
    if (comment) *comment = '\0'; // skip comments

    first = linebuf;
    last = linebuf + strlen(linebuf) - 1;
    while (*first == ' ' || *first == '\t') ++first;
    
    if (last > first) 
        while (*last == ' ' || *last == '\t' || *last == '\n') --last;

    tmp = realloc(*line,strlen(first)+1);
    if (tmp == NULL) {
        free(*line);
        *line = NULL;
        return mem_fail();
    }
    *line = tmp;
    
    if (first >= last) {
        **line = '\0';
        return 0;
    }
    
    *(last+1) = '\0';
    strcpy(*line,first);
    
    return strlen(*line);
}

/**
    Stores a line stripped of comments in the provided string. Returns
    the length of the line on success, FREESASA_FAIL if malloc/realloc
    errors.
 */
static int
next_line(char **line,
          FILE *fp) 
{
    char *linebuf = NULL;
    size_t len = 0;
    int ret;

    ret = getline(&linebuf,&len,fp);

    if (ret >= 0) ret = strip_line(line,linebuf);
    else ret = FREESASA_FAIL;

    free(linebuf);
    
    return ret;
}

int
freesasa_classifier_parse_class(const char *name)
{
#if HAVE_STRNCASECMP
    if (strncasecmp(name, "apolar", 6) == 0) {
        return FREESASA_ATOM_APOLAR;
    } else if (strncasecmp(name, "polar", 5) == 0) {
        return FREESASA_ATOM_POLAR;
    } else {
        return fail_msg("Only atom classes allowed are 'polar' and 'apolar'"
                        " (case insensitive)");
    }
#else
    if (strncmp(name, "apolar", 6) == 0) {
        return FREESASA_ATOM_APOLAR;
    } else if (strncmp(name, "polar", 5) == 0) {
        return FREESASA_ATOM_POLAR;
    } else {
        return fail_msg("Only atom classes allowed are 'polar' and 'apolar'");
    }
#endif
}

/**
    Add type. Returns the index of the new type on success,
    FREESASA_FAIL if realloc/strdup fails, FREESASA_WARN if type
    already known (ignore duplicates).
 */
int
freesasa_classifier_add_type(struct classifier_types *types,
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
        return freesasa_warn("Ignoring duplicate configuration entry for '%s'.", type_name);
    
    the_class = freesasa_classifier_parse_class(class_name);
    if (the_class == FREESASA_FAIL) return fail_msg("");
    
    if ((types->name = realloc(tn, sizeof(char*)*n)) == NULL) {
        types->name = tn;
        return mem_fail();
    }
    
    if ((types->type_radius = realloc(tr, sizeof(double)*n)) == NULL) {
        types->type_radius = tr;
        return mem_fail();
    }
    
    if ((types->type_class = realloc(tc, sizeof(int) * n)) == NULL) {
        types->type_class = tc;
        return mem_fail();
    }
    
    if ((types->name[n-1] = strdup(type_name)) == NULL) {
        return mem_fail();
    }
        
    types->n_types++;
    types->type_radius[types->n_types-1] = r;
    types->type_class[types->n_types-1] = the_class;
    return types->n_types-1;
}

/**
    Read a line specifying a type, store it in the config. Returns
    warning for duplicates, failures for syntax errors or memory
    allocation errors.
 */
static int
read_types_line(struct classifier_types *types,
                const char* line) 
{
    size_t blen=101;
    char buf1[blen], buf2[blen];
    int the_type;
    double r;
    if (sscanf(line,"%s %lf %s",buf1,&r,buf2) == 3) {
        the_type = freesasa_classifier_add_type(types, buf1, buf2, r);
        if (the_type == FREESASA_FAIL) return fail_msg("");
        if (the_type == FREESASA_WARN) return FREESASA_WARN;
    } else {
        return freesasa_fail("in %s(): Could not parse line '%s' in configuration, "
                             "expecting triplet of type 'TYPE [RADIUS] CLASS' for "
                             "example 'C_ALI 2.00 apolar'",
                             __func__, line);
    }
    return FREESASA_SUCCESS;
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
    char *line = NULL;
    int ret = FREESASA_SUCCESS, nl;
    size_t blen=101;
    char buf[blen];
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    if (fscanf(input,"%s",buf) == 0) return FREESASA_FAIL;
    assert(strcmp(buf,"types:") == 0);
    while (ftell(input) < fi.end) { 
        nl = next_line(&line,input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) {ret = nl; break; };
        ret = read_types_line(types,line);
        if (ret == FREESASA_FAIL) break;
    }
    free(line);
    return ret;
}

/**
    Add atom to residue. Returns index of the new atom on
    success. FREESASA_FAIL if memory allocation fails. FREESASA_WARN
    if the atom has already been added.
 */
int
freesasa_classifier_add_atom(struct classifier_residue *res,
                             const char *name,
                             double radius,
                             int the_class)
{
    int n;
    char **an = res->atom_name;
    double *ar = res->atom_radius;
    freesasa_atom_class *ac = res->atom_class;

    if (find_string(res->atom_name, name, res->n_atoms) >= 0)
        return freesasa_warn("Ignoring duplicate configuration entry for atom '%s %s'", 
                             res->name, name);
    n = res->n_atoms+1;

    if ((res->atom_name = realloc(res->atom_name,sizeof(char*)*n)) == NULL) {
        res->atom_name = an;
        return mem_fail();
    }
    if ((res->atom_radius = realloc(res->atom_radius,sizeof(double)*n)) == NULL) {
        res->atom_radius = ar;
        return mem_fail();
    }
    if ((res->atom_class = realloc(res->atom_class,sizeof(int)*n)) == NULL) {
        res->atom_class = ac;
        return mem_fail();
    }
    if ((res->atom_name[n-1] = strdup(name)) == NULL) 
        return mem_fail();

    ++res->n_atoms;
    res->atom_radius[n-1] = radius;
    res->atom_class[n-1] = the_class;

    return n-1;
}

/**
    Add residue to config. If the residue already exists, it returns
    the index of that residue, else it returns the index of the new
    residue. Returns FREESASA_FAILURE if realloc/strdup fails.
 */
int
freesasa_classifier_add_residue(struct freesasa_classifier *c,
                                const char* name)
{
    char **rn = c->residue_name;
    struct classifier_residue **cr = c->residue;
    int res = find_string(c->residue_name, name, c->n_residues);

    if (res >= 0) return res;

    res = c->n_residues + 1;
    if ((c->residue_name = realloc(rn, sizeof(char*) * res)) == NULL) {
        c->residue_name = rn;
        return mem_fail();
    }
    if ((c->residue = realloc(cr, sizeof(struct classifier_residue *) * res)) == NULL) {
        c->residue = cr;
        return mem_fail();
    }
    if ((c->residue[res-1] = freesasa_classifier_residue_new(name)) == NULL) {
        return mem_fail();
    }
    ++c->n_residues;
    c->residue_name[res-1] = c->residue[res-1]->name;
    return res-1;
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
                const char* line)
{
    size_t blen=100;
    char buf1[blen], buf2[blen], buf3[blen];
    int res, type, atom;
    if (sscanf(line,"%s %s %s",buf1,buf2,buf3) == 3) {
        type = find_string(types->name, buf3, types->n_types);
        if (type < 0) 
            return freesasa_fail("Unknown atom type '%s' in configuration, line '%s'",
                                 buf3, line);
        res = freesasa_classifier_add_residue(c, buf1);
        if (res == FREESASA_FAIL) return fail_msg("");
        atom = freesasa_classifier_add_atom(c->residue[res],
                                            buf2,
                                            types->type_radius[type],
                                            types->type_class[type]);

        if (atom == FREESASA_FAIL) return fail_msg("");
        if (atom == FREESASA_WARN) return FREESASA_WARN;
        
    } else {
        return freesasa_fail("in %s(): Could not parse configuration, line '%s', "
                             "expecting triplet of type "
                             "'RESIDUE ATOM CLASS', for example 'ALA CB C_ALI'.",
                             __func__, line);
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
    size_t blen=100;
    char *line = NULL, buf[blen];
    int ret = FREESASA_SUCCESS, nl;
    fseek(input, fi.begin, SEEK_SET);
    // read command (and discard)
    if (fscanf(input, "%s", buf) == 0) return FREESASA_FAIL;
    assert(strcmp(buf, "atoms:") == 0);
    while (ftell(input) < fi.end) { 
        nl = next_line(&line, input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) return fail_msg("");
        ret = read_atoms_line(c, types, line);
        if (ret == FREESASA_FAIL) break;
    }
    free(line);

    return ret;
}

static struct freesasa_classifier*
read_config(FILE *input) 
{
    assert(input);
    struct file_range types_section, atoms_section; 
    struct freesasa_classifier *classifier = NULL;
    struct classifier_types *types = NULL;

    if (!(types = freesasa_classifier_types_new()))
        goto cleanup;
    if (!(classifier = freesasa_classifier_new()))
        goto cleanup;
    if (check_file(input, &types_section, &atoms_section))
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
    *res = find_string(c->residue_name,"ANY",c->n_residues);
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
          int* res,
          int* atom)
{
    *atom = -1;
    *res = find_string(c->residue_name, res_name, c->n_residues);
    if (*res < 0) {
        find_any(c, atom_name, res, atom);
    } else {        
        const struct classifier_residue *residue = c->residue[*res];
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
    assert(classifier); assert(res_name); assert(atom_name);
    
    int res, atom, status;

    status = find_atom(classifier, res_name, atom_name, &res,&atom);
    if (status == FREESASA_SUCCESS)
        return classifier->residue[res]->atom_radius[atom];
    return -1.0;
}

freesasa_atom_class
freesasa_classifier_class(const freesasa_classifier *classifier,
                          const char *res_name, 
                          const char *atom_name)
{
    assert(classifier); assert(res_name); assert(atom_name);
    int res, atom, status;

    status = find_atom(classifier, res_name, atom_name, &res, &atom);
    if (status == FREESASA_SUCCESS)
        return classifier->residue[res]->atom_class[atom];
    return FREESASA_ATOM_UNKNOWN;
}

const char*
freesasa_classifier_class2str(const freesasa_classifier *classifier,
                              int the_class)
                              
{
    assert(classifier);
    switch (the_class) {
    case FREESASA_ATOM_APOLAR:
        return "Apolar";
    case FREESASA_ATOM_POLAR:
        return "Polar";
    case FREESASA_ATOM_UNKNOWN:
        return "Unknown";
    }
    return NULL;
}

freesasa_nodearea
freesasa_classifier_classify_result(const freesasa_classifier *classifier,
                                    const freesasa_structure *structure,
                                    const freesasa_result *result)
{
    freesasa_nodearea area = {"whole-structure", 0, 0, 0, 0, 0};
    if (classifier == NULL) classifier = &freesasa_default_classifier;
    freesasa_range_nodearea(&area, structure, result, classifier,
                           0, freesasa_structure_n(structure) - 1);
    return area;
}

static
freesasa_classifier*
classifier_from_file(FILE *file, const char *name)
{
    struct freesasa_classifier *classifier = read_config(file);
    if (classifier == NULL) {
        fail_msg("");
        return NULL;
    }
    classifier->name = strdup(name);
    if (classifier->name == NULL) {
        mem_fail();
        return NULL;
    }
    return classifier;
}

freesasa_classifier*
freesasa_classifier_from_file(FILE *file)
{
    assert(file);
    return classifier_from_file(file, "from-unknown-file");
}

freesasa_classifier*
freesasa_classifier_from_filename(const char *filename)
{
    FILE *file = fopen(filename, "r");
    freesasa_classifier *c = NULL;
    if (file) {
        c = classifier_from_file(file, filename);
        fclose(file);
        if (c == NULL) fail_msg("");
    } else {
        freesasa_fail("Error: could not open file '%s'; %s",
                      filename, strerror(errno));
    }       
    return c;
}

const freesasa_nodearea *
freesasa_classifier_residue_reference(const freesasa_classifier *classifier,
                                      const char *res_name)                                      
{
    int res = find_string(classifier->residue_name, res_name, classifier->n_residues);
    if (res < 0) return NULL;
    
    return &classifier->residue[res]->max_area;
}

const char*
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
    // elements that actually occur in the regular amino acids and nucleotides
    {" H", 1.10}, {" C", 1.70}, {" N", 1.55}, {" O", 1.52}, {" P", 1.80}, {" S", 1.80}, {"SE", 1.90}, 
    // some others, just because there were readily available values
    {" F", 1.47}, {"CL", 1.75}, {"BR", 1.83}, {" I", 1.98},
    {"LI", 1.81}, {"BE", 1.53}, {" B", 1.92}, 
    {"NA", 2.27}, {"MG", 1.74}, {"AL", 1.84}, {"SI", 2.10}, 
    {" K", 2.75}, {"CA", 2.31}, {"GA", 1.87}, {"GE", 2.11}, {"AS", 1.85}, 
    {"RB", 3.03}, {"SR", 2.49}, {"IN", 1.93}, {"SN", 2.17}, {"SB", 2.06}, {"TE", 2.06}, 
};

double
freesasa_guess_radius(const char* symbol)
{
    assert(symbol);
    int n_symbol = sizeof(symbol_radius)/sizeof(struct symbol_radius);
    for (int i = 0; i < n_symbol; ++i) {
        if (strcmp(symbol,symbol_radius[i].symbol) == 0)
            return symbol_radius[i].radius;
    }
    return -1.0;
}

//! The residue types that are returned by freesasa_classify_residue()
enum residue {
    //Regular amino acids
    ALA=0, ARG, ASN, ASP,
    CYS, GLN, GLU, GLY,
    HIS, ILE, LEU, LYS, 
    MET, PHE, PRO, SER,
    THR, TRP, TYR, VAL,
    //some non-standard ones
    CSE, SEC, PYL, PYH,
    ASX, GLX,
    //residue unknown
    RES_UNK,
    //capping N- and C-terminal groups (usually HETATM)
    ACE, NH2,
    //DNA
    DA, DC, DG, DT,
    DU, DI,
    //RNA (avoid one-letter enums)
    RA, RC, RG, RU, RI, RT,
    //generic nucleotide
    NN
};

// Residue types, make sure this always matches the corresponding enum.
static const char *residue_names[] = {
    //amino acids
    "ALA","ARG","ASN","ASP",
    "CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    // non-standard amino acids
    "CSE","SEC","PYL","PYH", // SEC and PYL are standard names, CSE and PYH are found in some early files
    "ASX","GLX",
    "UNK",
    // capping groups
    "ACE","NH2",
    //DNA
    "DA","DC","DG","DT","DU","DI",
    //RNA
    "A","C","G","U","I","T",
    //General nucleotide
    "N"
};

int
freesasa_classify_n_residue_types()
{
    return NN+1;
}

int
freesasa_classify_residue(const char *res_name)
{
    int len = strlen(res_name);
    char cpy[len+1];
    sscanf(res_name, "%s", cpy);
    for (int i = ALA; i < freesasa_classify_n_residue_types(); ++i) {
        if (strcmp(cpy,residue_names[i]) == 0) return i;
    }
    return RES_UNK;
}

const char *
freesasa_classify_residue_name(int residue_type)
{
    assert(residue_type >= 0 && residue_type <= NN);
    return residue_names[residue_type];
}

int
freesasa_atom_is_backbone(const char *atom_name)
{
    int n = strlen(atom_name);
    char name[n+1];
    name[0] = '\0';
    sscanf(atom_name, "%s", name); //trim whitespace

    if (strlen(name) == 0) return 0;
    if (strcmp(name, "CA") == 0 ||
        strcmp(name, "N") == 0 ||
        strcmp(name, "O") == 0 ||
        strcmp(name, "C") == 0 ||
        strcmp(name, "OXT") == 0)
        return 1;
    return 0;
}

#if USE_CHECK
#include <check.h>
#include <math.h>

START_TEST (test_classifier)
{
    const char *strarr[] = {"A","B","C"};
    const char *line[] = {"# Bla"," # Bla","Bla # Bla"," Bla # Bla","#Bla #Alb"};
    char *dummy_str = NULL;
    struct classifier_types *types = freesasa_classifier_types_new();
    struct classifier_residue *residue_cfg = freesasa_classifier_residue_new("ALA");
    struct freesasa_classifier *clf = freesasa_classifier_new();

    freesasa_set_verbosity(FREESASA_V_SILENT);

    ck_assert_int_eq(find_string((char**)strarr,"A",3),0);
    ck_assert_int_eq(find_string((char**)strarr,"B",3),1);
    ck_assert_int_eq(find_string((char**)strarr,"C",3),2);
    ck_assert_int_eq(find_string((char**)strarr,"D",3),-1);
    ck_assert_int_eq(find_string((char**)strarr," C ",3),2);
    ck_assert_int_eq(find_string((char**)strarr,"CC",3),-1);

    ck_assert_int_eq(strip_line(&dummy_str,line[0]),0);
    ck_assert_int_eq(strip_line(&dummy_str,line[1]),0);
    ck_assert_int_eq(strip_line(&dummy_str,line[2]),3);
    ck_assert_str_eq(dummy_str,"Bla");
    ck_assert_int_eq(strip_line(&dummy_str,line[3]),3);
    ck_assert_str_eq(dummy_str,"Bla");
    ck_assert_int_eq(strip_line(&dummy_str,line[4]),0);

    ck_assert_int_eq(freesasa_classifier_parse_class("A"), FREESASA_FAIL);
#if HAVE_STRNCASECMP
    ck_assert_int_eq(freesasa_classifier_parse_class("POLAR"), FREESASA_ATOM_POLAR);
    ck_assert_int_eq(freesasa_classifier_parse_class("APOLAR"), FREESASA_ATOM_APOLAR);
#endif
    ck_assert_int_eq(freesasa_classifier_parse_class("polar"), FREESASA_ATOM_POLAR);
    ck_assert_int_eq(freesasa_classifier_parse_class("apolar"), FREESASA_ATOM_APOLAR);

    ck_assert_int_eq(types->n_types, 0);
    ck_assert_int_eq(freesasa_classifier_add_type(types,"a","A",1.0),FREESASA_FAIL);
    ck_assert_int_eq(freesasa_classifier_add_type(types,"a","polar",1.0),0);
    ck_assert_int_eq(freesasa_classifier_add_type(types,"b","apolar",2.0),1);
    ck_assert_int_eq(freesasa_classifier_add_type(types,"b","polar",1.0),FREESASA_WARN);
    ck_assert_int_eq(freesasa_classifier_add_type(types,"c","apolar",3.0),2);
    ck_assert_int_eq(types->n_types,3);
    ck_assert_str_eq(types->name[0],"a");
    ck_assert_str_eq(types->name[1],"b");
    ck_assert_str_eq(types->name[2],"c");
    ck_assert(fabs(types->type_radius[0]-1.0) < 1e-10);
    ck_assert(fabs(types->type_radius[1]-2.0) < 1e-10);
    ck_assert(fabs(types->type_radius[2]-3.0) < 1e-10);

    freesasa_classifier_types_free(types);
    types = freesasa_classifier_types_new();

    ck_assert_int_eq(read_types_line(types,""),FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types,"a"),FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types,"a 1.0"),FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types,"a b C"),FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types,"a 1.0 C"),FREESASA_FAIL);
    ck_assert_int_eq(read_types_line(types,"a 1.0 apolar"),FREESASA_SUCCESS);
    ck_assert_int_eq(read_types_line(types,"b 2.0 polar"),FREESASA_SUCCESS);
    ck_assert_int_eq(types->n_types,2);
    ck_assert_str_eq(types->name[0],"a");
    ck_assert_str_eq(types->name[1],"b");
    ck_assert(fabs(types->type_radius[0]-1.0) < 1e-10);
    ck_assert(fabs(types->type_radius[1]-2.0) < 1e-10);

    ck_assert_int_eq(freesasa_classifier_add_atom(residue_cfg,"C",1.0,0),0);
    ck_assert_int_eq(freesasa_classifier_add_atom(residue_cfg,"CB",2.0,0),1);
    ck_assert_int_eq(freesasa_classifier_add_atom(residue_cfg,"CB",2.0,0),FREESASA_WARN);
    ck_assert_str_eq(residue_cfg->atom_name[0],"C");
    ck_assert_str_eq(residue_cfg->atom_name[1],"CB");
    ck_assert(fabs(residue_cfg->atom_radius[0]-1.0) < 1e-10);
    ck_assert(fabs(residue_cfg->atom_radius[1]-2.0) < 1e-10);
    freesasa_classifier_residue_free(residue_cfg);

    ck_assert_int_eq(freesasa_classifier_add_residue(clf,"A"),0);
    ck_assert_int_eq(freesasa_classifier_add_residue(clf,"B"),1);
    ck_assert_int_eq(freesasa_classifier_add_residue(clf,"B"),1);
    ck_assert_int_eq(clf->n_residues,2);
    ck_assert_str_eq(clf->residue_name[0],"A");
    ck_assert_str_eq(clf->residue_name[1],"B");
    ck_assert_str_eq(clf->residue[0]->name,"A");

    freesasa_classifier_free(clf);
    clf = freesasa_classifier_new();

    ck_assert_int_eq(read_atoms_line(clf,types,"A A"),FREESASA_FAIL);
    ck_assert_int_eq(read_atoms_line(clf,types,"A A bla"),FREESASA_FAIL);
    ck_assert_int_eq(read_atoms_line(clf,types,"ALA CA a"),FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf,types,"ALA CB b"),FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf,types,"ARG CA a"),FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf,types,"ARG CB b"),FREESASA_SUCCESS);
    ck_assert_int_eq(read_atoms_line(clf,types,"ARG CG b"),FREESASA_SUCCESS);
    ck_assert_int_eq(clf->n_residues,2);
    ck_assert_str_eq(clf->residue_name[0],"ALA");
    ck_assert_str_eq(clf->residue_name[1],"ARG");
    ck_assert_int_eq(clf->residue[0]->n_atoms,2);
    ck_assert_str_eq(clf->residue[0]->atom_name[0],"CA");
    ck_assert_str_eq(clf->residue[0]->atom_name[1],"CB");
    ck_assert(fabs(clf->residue[0]->atom_radius[0]-1.0) < 1e-5);
    ck_assert(fabs(clf->residue[0]->atom_radius[1]-2.0) < 1e-5);

    freesasa_classifier_free(clf);
    freesasa_classifier_types_free(types);

    freesasa_set_verbosity(FREESASA_V_NORMAL);

    free(dummy_str);
}
END_TEST

TCase *
test_classifier_static()
{
    TCase *tc = tcase_create("classifier.c static");
    tcase_add_test(tc, test_classifier);

    return tc;
}

#endif // USE_CHECK
