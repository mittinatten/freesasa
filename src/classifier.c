/*
  Copyright Simon Mitternacht 2013-2015.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "freesasa.h"
#include "util.h"

/**
    In this file the concept class refers to polar/apolar and type to
    aliphatic/aromatic/etc. See the example configurations in share/.
 */

static const char *default_type_input[] = {
    "C_ALI 2.00 apolar ",
    "C_ARO 1.75 apolar",
    "C_CAR 1.55 polar",
    "N 1.55 polar",
    "O 1.40 polar", // carbo- and hydroxyl oxygen have the same radius in OONS
    "S 2.00 polar",
    "SE 1.90 polar",
};

static const char *default_atom_input[] = {
    "ANY C   C_CAR",
    "ANY O   O",
    "ANY CA  C_ALI",
    "ANY N   N",
    "ANY CB  C_ALI",
    "ANY OXT O",

    "ARG CG C_ALI",
    "ARG CD C_ALI",
    "ARG NE N",
    "ARG CZ C_ALI",
    "ARG NH1 N",
    "ARG NH2 N",

    "ASN CG  C_CAR",
    "ASN OD1 O",
    "ASN ND2 N",

    "ASP CG  C_CAR",
    "ASP OD1 O",
    "ASP OD2 O",

    "CYS SG  S",

    "GLN CG  C_ALI",
    "GLN CD  C_CAR",
    "GLN OE1 O",
    "GLN NE2 N",

    "GLU CG  C_ALI",
    "GLU CD  C_CAR",
    "GLU OE1 O",
    "GLU OE2 O",

    "HIS CG  C_ARO",
    "HIS ND1 N",
    "HIS CD2 C_ARO",
    "HIS NE2 N",
    "HIS CE1 C_ARO",

    "ILE CG1 C_ALI",
    "ILE CG2 C_ALI",
    "ILE CD1 C_ALI",

    "LEU CG  C_ALI",
    "LEU CD1 C_ALI",
    "LEU CD2 C_ALI",

    "LYS CG  C_ALI",
    "LYS CD  C_ALI",
    "LYS CE  C_ALI",
    "LYS NZ  N",

    "MET CG  C_ALI",
    "MET SD  S",
    "MET CE  C_ALI",

    "PHE CG  C_ARO",
    "PHE CD1 C_ARO",
    "PHE CD2 C_ARO",
    "PHE CE1 C_ARO",
    "PHE CE2 C_ARO",
    "PHE CZ  C_ARO",

    "PRO CB  C_ARO",
    "PRO CG  C_ARO",
    "PRO CD  C_ARO",

    "SER OG  O",

    "THR OG1 O",
    "THR CG2 C_ALI",

    "TRP CG  C_ARO",
    "TRP CD1 C_ARO",
    "TRP CD2 C_ARO",
    "TRP NE1 N",
    "TRP CE2 C_ARO",
    "TRP CE3 C_ARO",
    "TRP CZ2 C_ARO",
    "TRP CZ3 C_ARO",
    "TRP CH2 C_ARO",

    "TYR CG  C_ARO",
    "TYR CD1 C_ARO",
    "TYR CD2 C_ARO",
    "TYR CE1 C_ARO",
    "TYR CE2 C_ARO",
    "TYR CZ  C_ARO",
    "TYR OH  O",

    "VAL CG1 C_ALI",
    "VAL CG2 C_ALI",

    "SEC SE SE",
    // add more here
};

// count references to default classifier
static unsigned int default_classifier_refcount = 0;
static freesasa_classifier *default_classifier = NULL;

/**
    Struct to store information about the types-section in a user-config.
 */
struct types {
    int n_classes; //!< number of classes
    int n_types; //!< number of types
    char **name; //!< names of types
    double *type_radius; //!< radius of type
    int *type_class; //!< class of each type
    char **class_name; //!< name of each type
};

static const struct types empty_types = {0, 0, NULL, NULL, NULL, NULL};

/**
     Configuration info for each residue type.
 */
struct residue_cfg {
    int n_atoms; //!< Number of atoms
    char *name; //!< Name of residue
    char **atom_name; //!< Names of atoms
    double *atom_radius; //!< Atomic radii
    int *atom_class; //!< Classes of atoms
};

static const struct residue_cfg empty_residue = {0,NULL,NULL,NULL,NULL};

/**
    Stores a user-configuration as extracted from a configuration
    file. No info about types, since those are only a tool used in
    assigment of radii and classes.
    
    An array of the names of residues is stored directly in the struct
    to facilitate searching for residues. The class_name array should
    be a clone of that found in types (can be done bye
    config_t_copy_classes()).
 */
struct config {
    int n_residues; //!< Number of residues
    int n_classes; //!< Number of classes
    char **residue_name; //!< Names of residues
    char **class_name; //!< Names of classes
    struct residue_cfg **residue;
};

static const struct config empty_config = {0, 0, NULL, NULL, NULL};

static struct types*
types_new()
{
    struct types *t = malloc(sizeof(struct types));
    if (t == NULL) {
        mem_fail();
        return NULL;
    }
    *t = empty_types;
    return t;
}

static void
types_free(struct types* t)
{
    if (t == NULL) return;
    free(t->type_radius);
    free(t->type_class);

    if (t->name)
        for (int i = 0; i < t->n_types; ++i)
            free(t->name[i]);
    free(t->name);

    if (t->class_name)
        for (int i = 0; i < t->n_classes; ++i)
            free(t->class_name[i]);
    free(t->class_name);

    free(t);
}

static struct residue_cfg*
residue_cfg_new(const char* name)
{
    assert(strlen(name) > 0);
    struct residue_cfg *res = malloc(sizeof(struct residue_cfg));
    if (res == NULL) {
        mem_fail();
        return NULL;
    }
    *res = empty_residue;
    res->name = strdup(name);
    if (res->name == NULL) {
        mem_fail();
        free(res);
        return NULL;
    }
    return res;
}

static void
residue_cfg_free(struct residue_cfg* res)
{
    if (res == NULL) return;
    free(res->name);

    if (res->atom_name)
        for (int i = 0; i < res->n_atoms; ++i)
            free(res->atom_name[i]);
    free(res->atom_name);

    free(res->atom_radius);
    free(res->atom_class);

    free(res);
}

static struct config* 
config_new()
{
    struct config *cfg = malloc(sizeof(struct config));
    if (cfg == NULL) {
        mem_fail();
        return NULL;
    }
    *cfg = empty_config;
    return cfg;
}

static void
config_free(void *p)
{
    if (p == NULL) return;
    struct config *c = p;

    if (c->class_name)
        for (int i = 0; i < c->n_classes; ++i)
            free(c->class_name[i]);
    free(c->class_name);

    if (c->residue)
        for (int i = 0; i < c->n_residues; ++i)
            residue_cfg_free(c->residue[i]);
    free(c->residue);
    free(c->residue_name);

    free(c);
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
        if (array[i] == NULL) {freesasa_fail("%d %d %s %s",i,array_size,key,key_trimmed);}
        assert(array[i]);
        if (strcmp(array[i],key_trimmed) == 0) return i;
    }
    return -1;
}

/**
    Checks that input file has the required fields and locates the
    'types' and 'atoms' sections. No syntax checking. Return
    FREESASA_SUCCESS if file seems ok, FREESASA_FILE if either/both of
    the sections are missing.
 */
static int 
check_file(FILE *input,
           struct file_interval *types, 
           struct file_interval *atoms)
{
    assert(input); assert(types); assert(atoms);
    long last_tell;
    const int blen = 200;
    char buf[blen];
    struct file_interval *last_interval = NULL;
        
    last_tell = ftell(input);
    types->begin = atoms->begin = -1;
    while (fscanf(input,"%s",buf) > 0) {
        if (strcmp(buf,"types:") == 0) {
            types->begin = last_tell;
            if (last_interval) last_interval->end = last_tell;
            last_interval = types;
        }
        if (strcmp(buf,"atoms:") == 0) {
            atoms->begin = last_tell;
            if (last_interval) last_interval->end = last_tell;
            last_interval = atoms;
        }
        last_tell = ftell(input);
    }
    if (last_interval != NULL) { 
        last_interval->end = last_tell;
    } 
    rewind(input);

    if ((types->begin == -1) || 
        (atoms->begin == -1)) {
        return freesasa_fail("in %s(): Input configuration lacks (at least) one of "
                             "the entries 'types:' or "
                             "'atoms:'.", __func__);
    }
    
    return FREESASA_SUCCESS;
}

/**
   Removes comments and strips leading and trailing
   whitespace. Returns the length of the stripped line on success,
   FREESASA_FAIL if malloc/realloc fails.
 */
int
strip_line(char **line,
           const char *input) 
{
    char *linebuf = malloc(strlen(input)+1),
        *comment, *first, *last;
    if (linebuf == NULL) return mem_fail();
    
    strcpy(linebuf,input);
    comment = strchr(linebuf,'#');
    if (comment) *comment = '\0'; // skip comments
    
    first = linebuf;
    last = linebuf + strlen(linebuf) - 1;
    while (*first == ' ' || *first == '\t') ++first;
    
    if (last > first) 
        while (*last == ' ' || *last == '\t' || *last == '\n') --last;
    *line = realloc(*line,strlen(first)+1);
    if (*line == NULL) return mem_fail();
    
    if (first >= last) {
        **line = '\0';
        free(linebuf);
        return 0;
    }
    
    *(last+1) = '\0';
    strcpy(*line,first);
    free(linebuf);
    
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
    
    getline(&linebuf,&len,fp);
    ret = strip_line(line,linebuf);
    free(linebuf);
    
    return ret;
}
/**
    Add class to type-registry. Returns the index of the new class on
    success, FREESASA_FAILURE if realloc/strdup fails.
 */
static int
add_class(struct types *types,
          const char *name)
{
    int the_class = find_string(types->class_name, name, types->n_classes);
    if (the_class < 0) {
        types->n_classes++;
        if (!(types->class_name = realloc(types->class_name, sizeof(char*)*types->n_classes))){
            types->n_classes--;
            return mem_fail();
        }
        if (!(types->class_name[types->n_classes-1] = strdup(name))) {
            types->n_classes--;
            return mem_fail();
        }
        the_class = types->n_classes - 1;
    }
    return the_class;
}

/**
    Add type. Returns the index of the new type on success,
    FREESASA_FAIL if realloc/strdup fails, FREESASA_WARN if type
    already known (ignore duplicates).
 */
static int
add_type(struct types *types,
         const char *type_name,
         const char *class_name, 
         double r)
{
    int the_class;
    if (find_string(types->name, type_name, types->n_types) >= 0)
        return freesasa_warn("Ignoring duplicate entry for '%s'.", type_name);
    the_class = add_class(types,class_name);
    if (the_class == FREESASA_FAIL) return fail_msg("");

    types->n_types++;
    types->name = realloc(types->name,sizeof(char*)*types->n_types);
    types->type_radius = realloc(types->type_radius,sizeof(double)*types->n_types);
    types->type_class = realloc(types->type_class,sizeof(int) * types->n_types);
    if (!types->name || !types->type_radius || !types->type_class) {
        types->n_types--;
        return mem_fail();
    }
    if (!(types->name[types->n_types-1] = strdup(type_name))) {
        types->n_types--;
        return mem_fail();
    }
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
read_types_line(struct types *types,
                const char* line) 
{
    size_t blen=101;
    char buf1[blen], buf2[blen];
    int the_type;
    double r;
    if (sscanf(line,"%s %lf %s",buf1,&r,buf2) == 3) {
        the_type = add_type(types, buf1, buf2, r);
        if (the_type == FREESASA_FAIL) return fail_msg("");
        if (the_type == FREESASA_WARN) return FREESASA_WARN;
    } else {
        return freesasa_fail("in %s(): Could not parse line '%s', expecting triplet of type "
                             "'TYPE [RADIUS] CLASS' for example 'C_ALI 2.00 apolar'",
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
read_types(struct types *types,
           FILE *input,
           struct file_interval fi)
{
    char *line = NULL;
    int ret = FREESASA_SUCCESS, nl;
    size_t blen=101;
    char buf[blen];
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf);
    assert(strcmp(buf,"types:") == 0);
    while (ftell(input) < fi.end) { 
        nl = next_line(&line,input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) return FREESASA_FAIL;
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
static int
add_atom(struct residue_cfg *res,
         const char *name,
         double radius,
         int the_class)
{
    int n;
    if (find_string(res->atom_name, name, res->n_atoms) >= 0)
        return freesasa_warn("in %s(): Ignoring duplicate entry for atom '%s %s'", 
                             __func__, res->name, name);
    n = ++res->n_atoms;
    res->atom_name = realloc(res->atom_name,sizeof(char*)*n);
    res->atom_radius = realloc(res->atom_radius,sizeof(double)*n);
    res->atom_class = realloc(res->atom_class,sizeof(int)*n);
    if (!res->atom_name || !res->atom_radius || !res->atom_class) {
        --res->n_atoms;
        return mem_fail();
    }
    if (!(res->atom_name[n-1] = strdup(name))) {
        --res->n_atoms;
        return mem_fail();
    }
    res->atom_radius[n-1] = radius;
    res->atom_class[n-1] = the_class;
    return n-1;
}

/**
    Add residue to config. If the residue already exists, it returns
    the index of that residue, else it returns the index of the new
    residue. Returns FREESASA_FAILURE if realloc/strdup fails.
 */
static int
add_residue(struct config *config,
            const char* name)
{
    int res = find_string(config->residue_name, name, config->n_residues);
    if (res >= 0) return res;
    res = ++config->n_residues;
    config->residue_name = realloc(config->residue_name, sizeof(char*) * res);
    config->residue = realloc(config->residue, sizeof(struct residue_cfg) * res);
    if (!config->residue_name || !config->residue) {
        --config->n_residues;
        return mem_fail();
    }
    if (!(config->residue[res-1] = residue_cfg_new(name))) {
        --config->n_residues;
        return mem_fail();
    }
    config->residue_name[res-1] = config->residue[res-1]->name;
    return res-1;
}

/**
    Read a line specifying an atom, store it in the config. Use
    supplied types to add assign radius and class. Returns
    FREESASA_WARN for duplicates. Returns FREESASA_FAIL for syntax
    errors or memory allocation errors. FREESASA_SUCCESS else.
 */
static int
read_atoms_line(struct config *config,
                const struct types *types,
                const char* line)
{
    size_t blen=100;
    char buf1[blen], buf2[blen], buf3[blen];
    int res, type, atom;
    if (sscanf(line,"%s %s %s",buf1,buf2,buf3) == 3) {
        type = find_string(types->name, buf3, types->n_types);
        if (type < 0) 
            return freesasa_fail("Unknown atom type '%s' in line '%s'",buf3,line);
        res = add_residue(config,buf1);
        if (res == FREESASA_FAIL) return fail_msg("");
        atom = add_atom(config->residue[res],
                        buf2,
                        types->type_radius[type],
                        types->type_class[type]);

        if (atom == FREESASA_FAIL) return fail_msg("");
        if (atom == FREESASA_WARN) return FREESASA_WARN;
        
    } else {
        return freesasa_fail("in %s(): Could not parse line '%s', expecting triplet of type "
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
read_atoms(struct config *config,
           struct types *types,
           FILE *input,
           struct file_interval fi)
{
    size_t blen=100;
    char *line = NULL, buf[blen];
    int ret = FREESASA_SUCCESS, nl;
    fseek(input, fi.begin, SEEK_SET);
    // read command (and discard)
    fscanf(input, "%s", buf);
    assert(strcmp(buf, "atoms:") == 0);
    while (ftell(input) < fi.end) { 
        nl = next_line(&line, input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) return fail_msg("");
        ret = read_atoms_line(config, types, line);
        if (ret == FREESASA_FAIL) break;
    }
    free(line);

    return ret;
}

static int
config_copy_classes(struct config *config,
                    const struct types *types) 
{
    char **names = malloc(sizeof(char*)*types->n_classes);
    if (names == NULL) return mem_fail();
    
    for (int i = 0; i < types->n_classes; ++i) {
        assert(types->class_name[i]);
        names[i] = strdup(types->class_name[i]);
        if (names[i] == NULL) return mem_fail();
    }
    config->n_classes = types->n_classes;
    config->class_name = names;
    return FREESASA_SUCCESS;
}

static struct config*
read_config(FILE *input) 
{
    assert(input);
    struct file_interval types_section, atoms_section; 
    struct config *config;
    struct types *types;
    
    if (!(types = types_new())) 
        return NULL;
    if (!(config = config_new()))
        return NULL;
    if (check_file(input, &types_section, &atoms_section) != FREESASA_SUCCESS)
        return NULL;
    
    if (read_types(types, input, types_section)         == FREESASA_FAIL ||
        read_atoms(config, types, input, atoms_section) == FREESASA_FAIL ||
        config_copy_classes(config, types)         == FREESASA_FAIL) {
        types_free(types);
        config_free(config);
        return NULL;
    }
    types_free(types);
    
    return config;
}

/**
    See if an atom_name has been defined for the residue ANY (writes
    indices to the provided pointers).
 */
static void 
find_any(const struct config *config,
         const char *atom_name,
         int *res, int *atom)
{
    *res = find_string(config->residue_name,"ANY",config->n_residues);
    if (*res >= 0) {
        *atom = find_string(config->residue[*res]->atom_name,
                            atom_name,
                            config->residue[*res]->n_atoms); 
    }
}
/**
    Find the residue and atom index of an atom in the supplied
    configuration. Prints error and returns FREESASA_FAIL if not
    found.
 */
static int 
find_atom(const struct config *config, 
          const char *res_name,
          const char *atom_name,
          int* res,
          int* atom)
{
    *atom = -1;
    *res = find_string(config->residue_name,res_name,config->n_residues);
    if (*res < 0) {
        find_any(config,atom_name,res,atom);
    } else {        
        const struct residue_cfg *residue = config->residue[*res];
        *atom = find_string(residue->atom_name,atom_name,residue->n_atoms);
        if (*atom < 0) {
            find_any(config,atom_name,res,atom);
        }
    }
    if (*atom < 0) {
        //return freesasa_warn("in %s(): Unknown residue '%s' and/or atom '%s'.",
        //                     __func__, res_name, atom_name);
        return FREESASA_WARN;
    }
    return FREESASA_SUCCESS;
}

/** To be linked to a Classifier struct */
static double
user_radius(const char *res_name,
            const char *atom_name,
            const freesasa_classifier *classifier)
{
    assert(classifier); assert(res_name); assert(atom_name);
    
    int res, atom, status;
    const struct config *config = classifier->config;
    
    status = find_atom(config,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return config->residue[res]->atom_radius[atom];
    //freesasa_warn("in %s(): couldn't find radius of atom '%s %s'.",
    //              __func__, res_name, atom_name);
    return -1.0;
}

/** To be linked to a Classifier struct */
static int
user_class(const char *res_name, 
           const char *atom_name,
           const freesasa_classifier *classifier)
{
    assert(classifier); assert(res_name); assert(atom_name);
    int res, atom, status;
    const struct config* config = classifier->config;
    status = find_atom(config,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return config->residue[res]->atom_class[atom];
    return FREESASA_WARN; //freesasa_warn("in %s(): couldn't find classification of atom '%s %s'.",
    //__func__, res_name, atom_name);
}

/** To be linked to a Classifier struct */
static const char*
user_class2str(int the_class,
               const freesasa_classifier *classifier)
{
    assert(classifier);
    const struct config* config = classifier->config;
    if (the_class < 0 && the_class >= config->n_classes) return NULL;
    return config->class_name[the_class];
}
static freesasa_classifier*
init_classifier(struct config *config)
{
   freesasa_classifier* c = malloc(sizeof(freesasa_classifier));
    if (c == NULL) {
        mem_fail();
        return NULL;
    }

    c->config = config;
    c->n_classes = config->n_classes;
    c->radius = user_radius;
    c->sasa_class = user_class;
    c->class2str = user_class2str;
    c->free_config = config_free;

    return c;
}

freesasa_classifier*
freesasa_classifier_from_file(FILE *file)
{
    assert(file);

    struct config *config = read_config(file);
    if (config == NULL) {
        fail_msg("");
        return NULL;
    }
    return init_classifier(config);
}

void
freesasa_classifier_free(freesasa_classifier *classifier)
{
    if (classifier != NULL) {
        if (classifier->free_config != NULL &&
            classifier->config != NULL) {
            classifier->free_config(classifier->config);
        }
        free(classifier);
    }
}


static struct config *
default_config()
{
    struct types *types;
    struct config *config;
    const int n_types = sizeof(default_type_input)/sizeof(char *);
    const int n_atoms = sizeof(default_atom_input)/sizeof(char *);

    types = types_new();
    if (types == NULL) return NULL;
    for (int i = 0; i < n_types; ++i) {
        if (read_types_line(types, default_type_input[i])
            == FREESASA_FAIL) {
            // this should never happen
            fail_msg("Error setting up types for default classifier");
            types_free(types);
            return NULL;
        }
    }

    config = config_new();
    if (config == NULL) return NULL;
    for (int i = 0; i < n_atoms; ++i) {
        if (read_atoms_line(config, types, default_atom_input[i])
            == FREESASA_FAIL) {
            // this should never happen
            fail_msg("Error setting up atoms for default classifier");
            config_free(config);
            types_free(types);
            return NULL;
        }
    }

    if (config_copy_classes(config, types) == FREESASA_FAIL) {
        config_free(config);
        config = NULL;
    }
    
    types_free(types);
    return config;
}

static freesasa_classifier *
classifier_default_new()
{
    struct config *config = default_config();
    if (config == NULL) {
        fail_msg("");
        return NULL;
    }
    
    return init_classifier(config);
}

const freesasa_classifier *
freesasa_classifier_default_acquire()
{
    if (default_classifier == NULL) {
        assert(default_classifier_refcount == 0);
        default_classifier = classifier_default_new();
        if (default_classifier == NULL) {
            fail_msg("Failed to load default classifier.");
            return NULL;
        }
    }
    ++default_classifier_refcount;
    return default_classifier;
}

void
freesasa_classifier_default_release() {
    if (default_classifier_refcount == 0) return;
    if (--default_classifier_refcount == 0) {
        freesasa_classifier_free(default_classifier);
        default_classifier = NULL;
    }
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
