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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "freesasa.h"
#include "util.h"

/**
    Struct to store user-configurations for classification.

    The file format is documented in the Public API.

    The types (aliphatic/aromatic/...) are stored here, but are only
    used as intermediaries for constructing the classification, and
    are not used later on.
 */
typedef struct user_config {
    char **classes; // names of area classes (polar/subpolar/etc)
    char **types; // names of atom types (aliphatic/aromatic/etc)
    char **residues; // names of residue types
    char ***atoms; // names of atom types per residue
    int n_classes; // number of classes
    int n_types; // number of atom types
    int n_residues; // number of residue types
    int *n_atoms; // number of atoms per residue type
    int **atom_class; // sasa-class of each atom
    int *type_class; //sasa-class of each atom type
    double *type_radius; // radius of each atom type
    double **atom_radius; // radius of each atom in each residue
} user_config;

static user_config* user_config_new()
{
    user_config *config = malloc(sizeof(user_config));
    if (config == NULL) { 
        mem_fail(); 
        return NULL;
    }
    config->classes = NULL;
    config->types = NULL;
    config->residues = NULL;
    config->atoms = NULL;
    config->n_classes = 0;
    config->n_types = 0;
    config->n_residues = 0;
    config->n_atoms = NULL;
    config->atom_class = NULL;
    config->type_class = NULL;
    config->type_radius = NULL;
    config->atom_radius = NULL;
    return config;
}
static void user_config_free(void *p)
{
    user_config *config;
    if (p) {
        config = (user_config*)p;
        if (config->classes) for (int i = 0; i < config->n_classes; ++i) free(config->classes[i]);
        if (config->types)   for (int i = 0; i < config->n_types; ++i)   free(config->types[i]);
        if (config->atoms) 
            for (int i = 0; i < config->n_residues; ++i) {
                if (config->n_atoms[i] && config->atoms[i]) 
                    for (int j = 0; j < config->n_atoms[i]; ++j) 
                        free(config->atoms[i][j]);
                free(config->atoms[i]);
            }
        if (config->residues)    for (int i = 0; i < config->n_residues; ++i) free(config->residues[i]);
        if (config->atom_class)  for (int i = 0; i < config->n_residues; ++i) free(config->atom_class[i]);
        if (config->atom_radius) for (int i = 0; i < config->n_residues; ++i) free(config->atom_radius[i]);
        free(config->types);
        free(config->classes);
        free(config->atoms);
        free(config->residues);
        free(config->n_atoms);
        free(config->atom_class);
        free(config->type_class);
        free(config->type_radius);
        free(config->atom_radius);
        free(config);
    }
}

// check if array of strings has a string that matches key
static int find_string(char **array, const char *key, int array_size)
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
    'types' and 'atoms' sections. No syntax checking.
 */
static int check_file(FILE *input,
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

// removes comments and strips leading and trailing whitespace
int strip_line(char **line, const char *input) {
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

static int next_line(char **line, FILE *fp) {
    char *linebuf = NULL;
    size_t len = 0;
    int ret;
    
    getline(&linebuf,&len,fp);
    ret = strip_line(line,linebuf);
    free(linebuf);
    
    return ret;
}

static int read_types_line(user_config *config,
                           const char* line) 
{
    size_t blen=101;
    char buf1[blen], buf2[blen];
    int areac;
    double r;
    if (sscanf(line,"%s %lf %s",buf1,&r,buf2) == 3) {
        if (find_string(config->types, buf1, config->n_types) >= 0) {
            return freesasa_warn("Ignoring duplicate entry for '%s'.", buf1);
        }
        areac = find_string(config->classes, buf2, config->n_classes);
        if (areac < 0) {
            config->n_classes++;
            if (!(config->classes = realloc(config->classes, sizeof(char*)*config->n_classes)))
                return mem_fail();
            if (!(config->classes[config->n_classes-1] = strdup(buf2))) 
                return mem_fail();
            areac = config->n_classes - 1;
        }
        --config->n_types;
        config->types = realloc(config->types,sizeof(char*)*config->n_types);
        config->type_radius = realloc(config->type_radius,sizeof(double)*config->n_types);
        config->type_class = realloc(config->type_class,sizeof(int) * config->n_types);
        if (!config->types || !config->type_radius || !config->type_class ||
            !(config->types[config->n_types-1] = strdup(buf1))) {
            --config->n_types;
            return mem_fail();
        }
        config->type_radius[config->n_types-1] = r;
        config->type_class[config->n_types-1] = areac;
            return mem_fail();
    } else {
        return freesasa_fail("%s: Could not parse line '%s', "
                         "expecting triplet of type "
                         "'TYPE [RADIUS] CLASS' for example "
                         "'C_ALI 2.00 apolar'",
                         __func__, line);
    }
    return FREESASA_SUCCESS;
}

/**
    Reads info about types from the user config. Assoicates each type
    with a class and a radius in the config struct..
*/
static int read_types(user_config *config,
                      FILE *input,
                      struct file_interval fi)
{
    char *line = NULL;
    int ret, nl = FREESASA_SUCCESS;
    size_t blen=101;
    char buf[blen];
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf);
    assert(strcmp(buf,"types:") == 0);
    while (ftell(input) < fi.end) { 
        nl = next_line(&line,input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) return mem_fail();
        ret = read_types_line(config,line);
        if (ret == FREESASA_FAIL) break;
    }
    free(line);
    return ret;
}

static int add_new_residue(user_config *config,
                           const char* residue)
{
    int res = ++config->n_residues;
    config->residues = realloc(config->residues, sizeof(char*) * res);
    config->n_atoms = realloc(config->n_atoms, sizeof(int) * res);
    config->atoms = realloc (config->atoms, sizeof(char**) * res);
    config->atom_class = realloc(config->atom_class, sizeof(int*) * res);
    config->atom_radius = realloc(config->atom_radius, sizeof(int*) * res);
    if (!config->residues || !config->n_atoms || !config->atoms ||
        !config->atom_class || !config->atom_radius ||
        !(config->residues[res-1] = strdup(residue))) {
        --config->n_residues;
        return mem_fail();
    }
    config->n_atoms[res-1] = 0;
    config->atoms[res-1] = NULL;
    config->atom_class[res-1] = NULL;
    config->atom_radius[res-1] = NULL;
    return res-1;
}

static int add_new_atom(user_config *config,
                        const char *atom,
                        int res,
                        int type)
{
    assert(config); assert(atom); assert(res>=0); assert(type>=0);
    int n = ++config->n_atoms[res];
    config->atoms[res] = realloc(config->atoms[res],sizeof(char*)*n);
    config->atom_class[res] = realloc(config->atom_class[res],sizeof(int)*n);
    config->atom_radius[res] = realloc(config->atom_radius[res],sizeof(double)*n);
    if (!config->atoms || !config->atom_class || !config->atom_radius ||
        !(config->atoms[res][n-1] = strdup(atom))) {
        --config->n_atoms[res];
        return mem_fail();
    }
    config->atom_class[res][n-1] = config->type_class[type];
    config->atom_radius[res][n-1] = config->type_radius[type];
    return FREESASA_SUCCESS;
}

static int read_atoms_line(user_config *config,
                           const char* line)
{
    size_t blen=100;
    char buf1[blen], buf2[blen], buf3[blen];
    if (sscanf(line,"%s %s %s",buf1,buf2,buf3) == 3) {
        int res = find_string(config->residues, buf1, config->n_residues);
        int type = find_string(config->types, buf3, config->n_types);
        if (type < 0) 
            return freesasa_fail("Unknown atom type '%s'",buf3);
        if (res < 0) 
            res = add_new_residue(config,buf1);
        if (res == FREESASA_FAIL)
            return mem_fail();
        if (find_string(config->atoms[res], buf2, config->n_atoms[res]) >= 0)
            freesasa_warn("Ignoring duplicate entry '%s %s %s'", buf1, buf2, buf3);
        if (add_new_atom(config,buf2,res,type) == FREESASA_FAIL) 
            return mem_fail();
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
static int read_atoms(user_config *config,
                      FILE *input,
                      struct file_interval fi)
{
    size_t blen=100;
    char *line = NULL, buf[blen];
    int res, type, n, ret, nl;
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf);
    assert(strcmp(buf,"atoms:") == 0);
    while (ftell(input) < fi.end) { 
        nl = next_line(&line,input);
        if (nl == 0) continue;
        if (nl == FREESASA_FAIL) return mem_fail();
        ret = read_atoms_line(config,line);
        if (ret == FREESASA_FAIL) break;
    }
    free(line);
    return FREESASA_SUCCESS;
}

static user_config* read_config(FILE *input) 
{
    assert(input);
    struct file_interval types, atoms; 
    int ret1, ret2;
    user_config *config;
    
    if (check_file(input,&types, &atoms) != FREESASA_SUCCESS) 
        return NULL;
    config = user_config_new();

    if (read_types(config, input, types) == FREESASA_FAIL ||
        read_atoms(config, input, atoms) == FREESASA_FAIL) {
        user_config_free(config);
        return NULL;
    }

    return config;
}

static void find_any(const user_config *config,
                     const char *atom_name,
                     int *res, int *atom)
{
    *res = find_string(config->residues,"ANY",config->n_residues);
    if (*res >= 0) {
        *atom = find_string(config->atoms[*res],atom_name,config->n_atoms[*res]); 
    }
}

static int find_atom(const user_config *config, 
                     const char *res_name,
                     const char *atom_name,
                     int* res,
                     int* atom)
{
    *atom = -1;
    *res = find_string(config->residues,res_name,config->n_residues);
    if (*res < 0) {
        find_any(config,atom_name,res,atom);
    } else {
        *atom = find_string(config->atoms[*res],atom_name,config->n_atoms[*res]);
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

static double user_radius(const char *res_name,
                          const char *atom_name,
                          const freesasa_classifier *classifier)
{
    assert(classifier); assert(res_name); assert(atom_name);
    
    int res, atom, status;
    const user_config *config = classifier->config;
    
    status = find_atom(config,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return config->atom_radius[res][atom];
    freesasa_fail("%s: couldn't find radius of atom '%s %s'.",
                  __func__, res_name, atom_name);
    return -1.0;
}

static int user_class(const char *res_name, 
                      const char *atom_name,
                      const freesasa_classifier *classifier)
{
    assert(classifier); assert(res_name); assert(atom_name);
    int res, atom, status;
    const user_config* config = classifier->config;
    status = find_atom(config,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return config->atom_class[res][atom];
    return freesasa_fail("%s: couldn't find classification of atom '%s %s'.",
                         __func__, res_name, atom_name);
}

static const char* user_class2str(int the_class,
                                  const freesasa_classifier *classifier)
{
    assert(classifier);
    const user_config* config = classifier->config;
    if (the_class < 0 && the_class >= config->n_classes) return NULL;
    return config->classes[the_class];
}

freesasa_classifier* freesasa_classifier_from_file(FILE *file)
{
    assert(file);
    user_config *config;
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

void freesasa_classifier_free(freesasa_classifier *classifier)
{
    if (classifier != NULL) {
        if (classifier->free_config != NULL &&
            classifier->config != NULL) {
            classifier->free_config(classifier->config);
        }
        free(classifier);
    }
}
