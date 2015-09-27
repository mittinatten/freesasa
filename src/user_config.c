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

/**
    Struct to store information about the types-section in a user-config.
 */
struct user_types {
    int n_classes; //!< number of classes
    int n_types; //!< number of types
    char **name; //!< names of types
    double *type_radius; //!< radius of type
    int *type_class; //!< class of each type
    char **class_name; //!< name of each type
};

static const struct user_types empty_types = {0, 0, NULL, NULL, NULL, NULL};

/**
     Configuration info for each residue type.
 */
struct user_residue {
    int n_atoms; //!< Number of atoms
    char *name; //!< Name of residue
    char **atom_name; //!< Names of atoms
    double *atom_radius; //!< Atomic radii
    int *atom_class; //!< Classe of atoms
};

static const struct user_residue empty_residue = {0,NULL,NULL,NULL,NULL};

/**
    Stores a user-configuration as extracted from a configuration
    file. No info about types, since those are only a tool used in
    assigment of radii and classes.
    
    An array of the names of residues is stored directly in the struct
    to facilitate searching for residues. The class_name array should
    be a clone of that found in user_types (can be done bye
    user_config_copy_classes()).
 */
struct user_config {
    int n_residues; //!< Number of residues
    int n_classes; //!< Number of classes
    char **residue_name; //!< Names of residues
    char **class_name; //!< Names of classes
    struct user_residue **residue;
};

static const struct user_config empty_config = {0, 0, NULL, NULL, NULL};

static struct user_types*
user_types_new()
{
    struct user_types *t = malloc(sizeof(struct user_types));
    if (t == NULL) { mem_fail(); return NULL; }
    *t = empty_types;
    return t;
}

static void
user_types_free(struct user_types* t)
{
    if (t == NULL) return;
    free(t->type_radius);
    free(t->type_class);
    if (t->name) for (int i = 0; i < t->n_types; ++i) free(t->name[i]);
    if (t->class_name) for (int i = 0; i < t->n_classes; ++i) free(t->class_name[i]);
    free(t->name);
    free(t->class_name);
    free(t);
}

static struct user_residue*
user_residue_new(const char* name)
{
    assert(strlen(name) > 0);
    struct user_residue *res = malloc(sizeof(struct user_residue));
    if (res == NULL) { mem_fail(); return NULL; }
    *res = empty_residue;
    res->name = strdup(name);
    return res;
}

static void
user_residue_free(struct user_residue* res)
{
    if (res == NULL) return;
    free(res->name);
    if (res->atom_name) for (int i = 0; i < res->n_atoms; ++i)  free(res->atom_name[i]);
    free(res->atom_name);
    free(res->atom_radius);
    free(res->atom_class);
    free(res);
}

static struct user_config* 
user_config_new()
{
    struct user_config *cfg = malloc(sizeof(struct user_config));
    if (cfg == NULL) { mem_fail(); return NULL; }
    *cfg = empty_config;
    return cfg;
}

static void
user_config_free(void *p)
{
    if (p == NULL) return;
    struct user_config *c = p;
    if (c->class_name) for (int i = 0; i < c->n_classes; ++i) free(c->class_name[i]);
    if (c->residue) for (int i = 0; i < c->n_residues; ++i) user_residue_free(c->residue[i]);
    free(c->residue_name);
    free(c->class_name);
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
    char key_trimmed[n];

    // remove trailing and leading whitespace
    sscanf(key,"%s",key_trimmed);
    for (int i = 0; i < array_size; ++i) {
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
        return freesasa_fail("%s: Input configuration lacks (at least) one of "
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
add_class(struct user_types *types,
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
    FREESASA_FAILURE if realloc/strdup fails, FREESASA_WARN if type
    already known (ignore duplicates).
 */
static int
add_type(struct user_types *types,
         const char *type_name,
         const char *class_name, 
         double r)
{
    int the_class;
    if (find_string(types->name, type_name, types->n_types) >= 0)
        return freesasa_warn("Ignoring duplicate entry for '%s'.", type_name);
    the_class = add_class(types,class_name);
    if (the_class == FREESASA_FAIL) return freesasa_fail(__func__);

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
read_types_line(struct user_types *types,
                const char* line) 
{
    size_t blen=101;
    char buf1[blen], buf2[blen];
    int the_type;
    double r;
    if (sscanf(line,"%s %lf %s",buf1,&r,buf2) == 3) {
        the_type = add_type(types, buf1, buf2, r);
        if (the_type == FREESASA_FAIL) return freesasa_fail(__func__);
        if (the_type == FREESASA_WARN) return FREESASA_WARN;
    } else {
        return freesasa_fail("%s: Could not parse line '%s', expecting triplet of type "
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
read_types(struct user_types *types,
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
add_atom(struct user_residue *res,
         const char *name,
         double radius,
         int the_class)
{
    int n;
    if (find_string(res->atom_name, name, res->n_atoms) >= 0)
        return freesasa_warn("%s: Ignoring duplicate entry for atom '%s %s'", 
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
add_residue(struct user_config *config,
            const char* name)
{
    int res = find_string(config->residue_name, name, config->n_residues);
    if (res >= 0) return res;
    res = ++config->n_residues;
    config->residue_name = realloc(config->residue_name, sizeof(char*) * res);
    config->residue = realloc(config->residue, sizeof(struct user_residue) * res);
    if (!config->residue_name || !config->residue) {
        --config->n_residues;
        return mem_fail();
    }
    if (!(config->residue[res-1] = user_residue_new(name))) {
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
read_atoms_line(struct user_config *config,
                const struct user_types *types,
                const char* line)
{
    size_t blen=100;
    char buf1[blen], buf2[blen], buf3[blen];
    int res, type, atom;
    if (sscanf(line,"%s %s %s",buf1,buf2,buf3) == 3) {
        type = find_string(types->name, buf3, types->n_types);
        if (type < 0) return freesasa_fail("Unknown atom type '%s' in line '%s'",buf3,line);
        res = add_residue(config,buf1);
        if (res == FREESASA_FAIL) return freesasa_fail(__func__);
        atom = add_atom(config->residue[res], buf2, types->type_radius[type], types->type_class[type]);
        if (atom == FREESASA_FAIL) return freesasa_fail(__func__);
        if (atom == FREESASA_WARN) return FREESASA_WARN;
        
    } else {
        return freesasa_fail("%s: Could not parse line '%s', expecting triplet of type "
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
read_atoms(struct user_config *config,
           struct user_types *types,
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
        if (nl == FREESASA_FAIL) return freesasa_fail(__func__);
        ret = read_atoms_line(config, types, line);
        if (ret == FREESASA_FAIL) break;
    }
    free(line);

    return ret;
}

static int
user_config_copy_classes(struct user_config *config,
                         const struct user_types *types) 
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

static struct user_config*
read_config(FILE *input) 
{
    assert(input);
    struct file_interval types_section, atoms_section; 
    struct user_config *config;
    struct user_types *types;
    
    if (!(types = user_types_new())) 
        return NULL;
    if (!(config = user_config_new()))
        return NULL;
    if (check_file(input, &types_section, &atoms_section) != FREESASA_SUCCESS)
        return NULL;
    
    if (read_types(types, input, types_section)         == FREESASA_FAIL ||
        read_atoms(config, types, input, atoms_section) == FREESASA_FAIL ||
        user_config_copy_classes(config, types)         == FREESASA_FAIL) {
        user_types_free(types);
        user_config_free(config);
        return NULL;
    }
    user_types_free(types);
    
    return config;
}

/**
    See if an atom_name has been defined for the residue ANY (writes
    indices to the provided pointers).
 */
static void 
find_any(const struct user_config *config,
         const char *atom_name,
         int *res, int *atom)
{
    *res = find_string(config->residue_name,"ANY",config->n_residues);
    if (*res >= 0) {
        *atom = find_string(config->residue[*res]->atom_name,atom_name,config->residue[*res]->n_atoms); 
    }
}
/**
    Find the residue and atom index of an atom in the supplied
    configuration. Prints error and returns FREESASA_FAIL if not
    found.
 */
static int 
find_atom(const struct user_config *config, 
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
        const struct user_residue *residue = config->residue[*res];
        *atom = find_string(residue->atom_name,atom_name,residue->n_atoms);
        if (*atom < 0) {
            find_any(config,atom_name,res,atom);
        }
    }
    if (*atom < 0) {
        return freesasa_fail("%s: Unknown residue '%s' and/or atom '%s'.",
                             __func__, res_name, atom_name);
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
    const struct user_config *config = classifier->config;
    
    status = find_atom(config,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return config->residue[res]->atom_radius[atom];
    freesasa_fail("%s: couldn't find radius of atom '%s %s'.",
                  __func__, res_name, atom_name);
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
    const struct user_config* config = classifier->config;
    status = find_atom(config,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return config->residue[res]->atom_class[atom];
    return freesasa_fail("%s: couldn't find classification of atom '%s %s'.",
                         __func__, res_name, atom_name);
}

/** To be linked to a Classifier struct */
static const char*
user_class2str(int the_class,
               const freesasa_classifier *classifier)
{
    assert(classifier);
    const struct user_config* config = classifier->config;
    if (the_class < 0 && the_class >= config->n_classes) return NULL;
    return config->class_name[the_class];
}

freesasa_classifier*
freesasa_classifier_from_file(FILE *file)
{
    assert(file);
    struct user_config *config;
    freesasa_classifier* c = malloc(sizeof(freesasa_classifier));
    if (c == NULL) { mem_fail(); return NULL; }
    config = read_config(file);
    if (config == NULL) return NULL;

    c->radius = user_radius;
    c->sasa_class = user_class;
    c->class2str = user_class2str;
    c->free_config = user_config_free;
    c->config = config;
    c->n_classes = config->n_classes;
    return c;
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
