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

extern int freesasa_fail(const char *format,...);
extern int freesasa_warn(const char *format,...);

struct file_interval {long begin; long end;};

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
    int *atom_type_class; //sasa-class of each atom type
    double *atom_type_radius; // radius of each atom type
    double **atom_radius; // radius of each atom in each residue
} user_config;
///////////////////////////
// User-provided classes //
///////////////////////////


static user_config* user_config_new()
{
    user_config *config = malloc(sizeof(user_config));
    config->classes = NULL;
    config->types = NULL;
    config->residues = NULL;
    config->atoms = NULL;
    config->n_classes = 0;
    config->n_types = 0;
    config->n_residues = 0;
    config->n_atoms = NULL;
    config->atom_class = NULL;
    config->atom_type_class = NULL;
    config->atom_type_radius = NULL;
    config->atom_radius = NULL;
    return config;
}

static int find_string(char **array, const char *key, int array_size)
{
    assert(key);
    if (array == NULL || array_size == 0) return -1;

    int n = strlen(key);
    char key_trimmed[n];

    sscanf(key,"%s",key_trimmed);
    for (int i = 0; i < array_size; ++i) {
        assert(array[i]);
        if (strcmp(array[i],key_trimmed) == 0) return i;
    }
    return -1;
} 

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
    char *linebuf = malloc(strlen(input)), 
        *comment, *first, *last;
    
    strcpy(linebuf,input);
    comment = strchr(linebuf,'#');
    if (comment) *comment = '\0'; // skip comments
    
    first = linebuf;
    last = linebuf + strlen(linebuf) - 1;
    while (*first == ' ' || *first == '\t') ++first;
    
    if (last > first) 
        while (*last == ' ' || *last == '\t' || *last == '\n') --last;
    *line = realloc(*line,strlen(first));
    
    if (first >= last) {
        **line = '\0';
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

static int read_types(user_config *config,
                      FILE *input,
                      struct file_interval fi)
{
    size_t blen=100;
    char *line = NULL, buf1[blen], buf2[blen];
    double r;
    int areac;
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf1);
    assert(strcmp(buf1,"types:") == 0);
    while (ftell(input) < fi.end) { 
        if (next_line(&line,input) == 0) continue;
        if (sscanf(line,"%s %lf %s",buf1,&r,buf2) == 3) {
            if (find_string(config->types, buf1, config->n_types) >= 0) {
                freesasa_warn("Ignoring duplicate entry for '%s'.", buf1);
                continue;
            }
            areac = find_string(config->classes, buf2, config->n_classes);
            if (areac < 0) {
                config->n_classes++;
                config->classes = realloc(config->classes,
                                                sizeof(char*) * config->n_classes);
                config->classes[config->n_classes-1] = strdup(buf2);
                areac = config->n_classes - 1;
            }
            config->n_types++;
            config->types = realloc(config->types,
                                            sizeof(char*) * config->n_types);
            config->types[config->n_types-1] = strdup(buf1);
            config->atom_type_radius = realloc(config->atom_type_radius,
                                                 sizeof(double) * config->n_types);
            config->atom_type_radius[config->n_types-1] = r;
            config->atom_type_class = realloc(config->atom_type_class,
                                                 sizeof(int) * config->n_types);
            config->atom_type_class[config->n_types-1] = areac;
        } else {
            return freesasa_fail("%s: Could not parse line '%s', "
                                 "expecting triplet of type "
                                 "'TYPE [RADIUS] CLASS' for example "
                                 "'C_ALI 2.00 apolar'.",
                                 __func__, line);
        }
    }
    free(line);
    return FREESASA_SUCCESS;
}

static int read_atoms(user_config *config,
                      FILE *input,
                      struct file_interval fi)
{
    size_t blen=100;
    char *line = NULL, buf1[blen], buf2[blen], buf3[blen];
    double r;
    int res, type, n;
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf1);
    assert(strcmp(buf1,"atoms:") == 0);
    while (ftell(input) < fi.end) { 
        if (next_line(&line,input) == 0) continue;
        if (sscanf(line,"%s %s %s",buf1,buf2,buf3) == 3) {
            res = find_string(config->residues, buf1, config->n_residues);
            type = find_string(config->types, buf3, config->n_types);
            if (type < 0) {
                free(line);
                return freesasa_fail("Unknown atom type '%s'",buf3);
            }
            if (res < 0) {
                config->n_residues++;
                res = config->n_residues - 1;
                config->residues = realloc(config->residues,
                                            sizeof(char*) * config->n_residues);
                config->n_atoms = realloc(config->n_atoms,
                                           sizeof(int) * config->n_residues);
                config->atoms = realloc (config->atoms,
                                          sizeof(char**) * config->n_residues);
                config->atom_class = realloc(config->atom_class,
                                              sizeof(int*) * config->n_residues);
                config->atom_radius = realloc(config->atom_radius,
                                               sizeof(int*) * config->n_residues);
                config->residues[res] = strdup(buf1);
                config->n_atoms[res] = 0;
                config->atoms[res] = NULL;
                config->atom_class[res] = NULL;
                config->atom_radius[res] = NULL;
            } 
            if (find_string(config->atoms[res],buf2,config->n_atoms[res]) >= 0) {
                freesasa_warn("Ignoring duplicate entry '%s %s %s'", buf1, buf2, buf3);
                continue;
            }
            fflush(stdout);
            n = ++config->n_atoms[res];
            config->atoms[res] = realloc(config->atoms[res],sizeof(char*)*n);
            config->atom_class[res] = realloc(config->atom_class[res],sizeof(int)*n);
            config->atom_radius[res] = realloc(config->atom_radius[res],sizeof(double)*n);
            config->atoms[res][n-1] = strdup(buf2);
            config->atom_class[res][n-1] = config->atom_type_class[type];
            config->atom_radius[res][n-1] = config->atom_type_radius[type];
        } else {
            return freesasa_fail("%s: Could not parse line '%s', "
                                 "expecting triplet of type "
                                 "'RESIDUE ATOM CLASS' for example "
                                 "'ALA CB C_ALI'.",
                                 __func__, line);
        }
    }
    free(line);
    return FREESASA_SUCCESS;
}
user_config* read_config(FILE *input) 
{
    assert(input);
    struct file_interval types, atoms; 
    int ret1, ret2;
    user_config *config;
    
    if (check_file(input,&types, &atoms) != FREESASA_SUCCESS) 
        return NULL;
    config = user_config_new();
    ret1 = read_types(config, input, types);
    ret2 = read_atoms(config, input, atoms);
    if (ret1 != FREESASA_SUCCESS || ret2 != FREESASA_SUCCESS) 
        return NULL;
    return config;
}

void user_config_free(void *p)
{
    user_config *config;
    if (p) {
        config = p;
        free(config->classes);
        free(config->types);
        for (int i = 0; i < config->n_residues; ++i) {
            free(config->residues[i]);
            for (int j = 0; j < config->n_atoms[i]; ++j) { 
                free(config->atoms[i][j]);
            }
            free(config->atoms[i]);
            free(config->atom_class[i]);
            free(config->atom_radius[i]);
        }
        free(config->residues);
        free(config->atoms);
        free(config->n_atoms);
        free(config->atom_class);
        free(config->atom_type_class);
        free(config->atom_type_radius);
        free(config->atom_radius);
        free(config);
    }
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
    
    freesasa_classifier* c = malloc(sizeof(freesasa_classifier));
    user_config* config = read_config(file);

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
    }
}
