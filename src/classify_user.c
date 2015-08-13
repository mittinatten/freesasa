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
#include "classify.h"

extern int freesasa_fail(const char *format,...);
extern int freesasa_warn(const char *format,...);

struct file_interval {long begin; long end;};

struct freesasa_classify {
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
};

extern int freesasa_trim_whitespace(char *target, const char *src,
                                    int length);

static freesasa_classify* freesasa_classify_new()
{
    freesasa_classify *classes = malloc(sizeof(freesasa_classify));
    classes->classes = NULL;
    classes->types = NULL;
    classes->residues = NULL;
    classes->atoms = NULL;
    classes->n_classes = 0;
    classes->n_types = 0;
    classes->n_residues = 0;
    classes->n_atoms = NULL;
    classes->atom_class = NULL;
    classes->atom_type_class = NULL;
    classes->atom_type_radius = NULL;
    classes->atom_radius = NULL;
    return classes;
}

static int find_string(char **array, const char *key, int array_size)
{
    assert(key);
    if (array == NULL || array_size == 0) return -1;
    int n = strlen(key);
    char key_trimmed[n];
    freesasa_trim_whitespace(key_trimmed, key, n);
    for (int i = 0; i < array_size; ++i) {
        assert(array[i]);
        if (strcmp(array[i],key_trimmed) == 0) return i;
    }
    return -1;
} 

static int freesasa_classify_check_file(FILE *input,
                                        struct file_interval *types, 
                                        struct file_interval *atoms)
{
    assert(input); assert(types); assert(atoms);
    long last_tell = ftell(input);
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
    last_interval->end = last_tell;
    rewind(input);
    
    if ((types->begin == -1) || 
        (atoms->begin == -1)) {
        return freesasa_fail("Input configuration lacks (at least) one of "
                             "the entries 'types:' or "
                             "'atoms:'.");
    }
    
    return FREESASA_SUCCESS;
}

static int freesasa_classify_read_types(freesasa_classify *classes,
                                               FILE *input,
                                               struct file_interval fi)
{
    size_t blen=100;
    char buf1[blen], buf2[blen];
    double r;
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf1);
    assert(strcmp(buf1,"types:") == 0);
    while (ftell(input) < fi.end) { 
        if (fscanf(input,"%s %lf %s",buf1,&r,buf2) > 0) {
            if (find_string(classes->types, buf1, classes->n_types) >= 0) {
                freesasa_warn("Ignoring duplicate entry '%s'.", buf1);
                continue;
            }
            int areac = find_string(classes->classes, buf2, classes->n_classes);
            if (areac < 0) {
                classes->n_classes++;
                classes->classes = realloc(classes->classes,
                                                sizeof(char*) * classes->n_classes);
                classes->classes[classes->n_classes-1] = strdup(buf2);
                areac = classes->n_classes - 1;
            }
            classes->n_types++;
            classes->types = realloc(classes->types,
                                            sizeof(char*) * classes->n_types);
            classes->types[classes->n_types-1] = strdup(buf1);
            classes->atom_type_radius = realloc(classes->atom_type_radius,
                                                 sizeof(double) * classes->n_types);
            classes->atom_type_radius[classes->n_types-1] = r;
            classes->atom_type_class = realloc(classes->atom_type_class,
                                                 sizeof(int) * classes->n_types);
            classes->atom_type_class[classes->n_types-1] = areac;
        }
    }
    return FREESASA_SUCCESS;
}

static int freesasa_classify_read_atoms(freesasa_classify *classes,
                                        FILE *input,
                                        struct file_interval fi)
{
    size_t blen=100;
    char buf1[blen], buf2[blen], buf3[blen];
    double r;
    fseek(input,fi.begin,SEEK_SET);
    // read command (and discard)
    fscanf(input,"%s",buf1);
    assert(strcmp(buf1,"atoms:") == 0);
    while (ftell(input) < fi.end) { 
        if (fscanf(input,"%s %s %s",buf1,buf2,buf3) > 0) {
            int res = find_string(classes->residues, buf1, classes->n_residues);
            int type = find_string(classes->types, buf3, classes->n_types);
            if (type < 0) {
                return freesasa_fail("Unknown atom type '%s'",buf3);
            }
            if (res < 0) {
                classes->n_residues++;
                res = classes->n_residues - 1;
                classes->residues = realloc(classes->residues,
                                            sizeof(char*) * classes->n_residues);
                classes->n_atoms = realloc(classes->n_atoms,
                                           sizeof(int) * classes->n_residues);
                classes->atoms = realloc (classes->atoms,
                                          sizeof(char**) * classes->n_residues);
                classes->atom_class = realloc(classes->atom_class,
                                              sizeof(int*) * classes->n_residues);
                classes->atom_radius = realloc(classes->atom_radius,
                                               sizeof(int*) * classes->n_residues);
                classes->residues[res] = strdup(buf1);
                classes->n_atoms[res] = 0;
                classes->atoms[res] = NULL;
                classes->atom_class[res] = NULL;
                classes->atom_radius[res] = NULL;
            } 
            if (find_string(classes->atoms[res],buf2,classes->n_atoms[res]) >= 0) {
                freesasa_warn("Ignoring duplicate entry '%s %s %s'", buf1, buf2, buf3);
                continue;
            }
            fflush(stdout);
            int n = ++classes->n_atoms[res];
            classes->atoms[res] = realloc(classes->atoms[res],sizeof(char*)*n);
            classes->atom_class[res] = realloc(classes->atom_class[res],sizeof(int)*n);
            classes->atom_radius[res] = realloc(classes->atom_radius[res],sizeof(double)*n);
            classes->atoms[res][n-1] = strdup(buf2);
            classes->atom_class[res][n-1] = classes->atom_type_class[type];
            classes->atom_radius[res][n-1] = classes->atom_type_radius[type];
        }
    }
    
    return FREESASA_SUCCESS;
}
freesasa_classify* freesasa_classify_user(FILE *input) 
{
    assert(input);
    struct file_interval types, atoms; 
    int result = freesasa_classify_check_file(input,&types, &atoms);
    if (result != FREESASA_SUCCESS) return NULL;
    freesasa_classify *classes = freesasa_classify_new();
    freesasa_classify_read_types(classes, input, types);
    freesasa_classify_read_atoms(classes, input, atoms);
    return classes;
}
freesasa_classify* freesasa_classify_user_clone(const freesasa_classify* source)
{
    freesasa_classify *copy =freesasa_classify_new();
    assert(copy);
    
    copy->n_classes = source->n_classes;
    copy->n_types = source->n_types;
    copy->n_residues = source->n_residues;

    copy->classes = malloc(sizeof(char*)*source->n_classes);
    copy->types = malloc(sizeof(char*)*source->n_types);
    copy->residues = malloc(sizeof(char*)*source->n_residues);
    copy->atoms = malloc(sizeof(char**)*source->n_residues);
    copy->n_atoms = malloc(sizeof(int)*source->n_residues);
    copy->atom_class = malloc(sizeof(int*)*source->n_residues);
    copy->atom_type_class = malloc(sizeof(int)*source->n_types);
    copy->atom_type_radius = malloc(sizeof(double)*source->n_types);
    copy->atom_radius = malloc(sizeof(double*)*source->n_residues);
    
    assert(copy->classes); assert(copy->types);
    assert(copy->residues); assert(copy->atoms);
    assert(copy->n_atoms); assert(copy->atom_class);
    assert(copy->atom_type_class); assert(copy->atom_type_radius);
    assert(copy->atom_radius);
    
    memcpy(copy->n_atoms,source->n_atoms,sizeof(int)*source->n_residues);
    memcpy(copy->atom_type_class,source->atom_type_class,sizeof(int)*source->n_types);
    memcpy(copy->atom_type_radius,source->atom_type_radius,sizeof(double)*source->n_types);

    for (int i = 0; i < source->n_classes; ++i)
        copy->classes[i] = strdup(source->classes[i]);
    for (int i = 0; i < source->n_types; ++i)
        copy->types[i] = strdup(source->types[i]);
    for (int i = 0; i < source->n_residues; ++i) {
        copy->residues[i] = strdup(source->residues[i]);
        copy->atoms[i] = malloc(sizeof(char*)*source->n_atoms[i]);
        copy->atom_class[i] = malloc(sizeof(int)*source->n_atoms[i]);
        copy->atom_radius[i] = malloc(sizeof(double)*source->n_atoms[i]);
        memcpy(copy->atom_class[i],source->atom_class[i],sizeof(int)*source->n_atoms[i]);
        memcpy(copy->atom_radius[i],source->atom_radius[i],sizeof(double)*source->n_atoms[i]);
        for (int j = 0; j < source->n_atoms[i]; ++j) {
            copy->atoms[i][j] = strdup(source->atoms[i][j]);
        }
    }
    return copy;
}

void freesasa_classify_user_free(freesasa_classify* classes)
{
    if (classes) {
        free(classes->classes);
        free(classes->types);
        for (int i = 0; i < classes->n_residues; ++i) {
            free(classes->residues[i]);
            for (int j = 0; j < classes->n_atoms[i]; ++j) { 
                free(classes->atoms[i][j]);
            }
            free(classes->atoms[i]);
            free(classes->atom_class[i]);
            free(classes->atom_radius[i]);
        }
        free(classes->residues);
        free(classes->atoms);
        free(classes->n_atoms);
        free(classes->atom_class);
        free(classes->atom_type_class);
        free(classes->atom_type_radius);
        free(classes->atom_radius);
    }
}

static void freesasa_classify_find_any(const freesasa_classify *classes,
                                      const char *atom_name,
                                      int *res, int *atom)
{
    *res = find_string(classes->residues,"ANY",classes->n_residues);
    if (*res >= 0) {
        *atom = find_string(classes->atoms[*res],atom_name,classes->n_atoms[*res]); 
    }
}

static int freesasa_classify_find_atom(const freesasa_classify *classes, 
                                       const char *res_name, const char *atom_name,
                                       int* res, int* atom)
{
    *atom = -1;
    *res = find_string(classes->residues,res_name,classes->n_residues);
    if (*res < 0) {
        freesasa_classify_find_any(classes,atom_name,res,atom);
    } else {
        *atom = find_string(classes->atoms[*res],atom_name,classes->n_atoms[*res]);
        if (*atom < 0) {
            freesasa_classify_find_any(classes,atom_name,res,atom);
        }
    }
    if (*atom < 0) {
        return freesasa_fail("%s: Unknown residue '%s' and/or atom '%s'.",
                             __func__, res_name, atom_name);
    }
    return FREESASA_SUCCESS;
}

double freesasa_classify_user_radius(const freesasa_classify *classes, 
                                     const char *res_name, const char *atom_name)
{
    assert(classes); assert(res_name); assert(atom_name);
    int res, atom, status;
    status = freesasa_classify_find_atom(classes,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return classes->atom_radius[res][atom];
    freesasa_fail("%s: couldn't find radius of atom '%s %s'.",
                  __func__, res_name, atom_name);
    return -1.0;
}
int freesasa_classify_user_n_classes(const freesasa_classify* classes)
{
    assert(classes);
    return classes->n_classes;
}

int freesasa_classify_user_class(const freesasa_classify *classes,
                                 const char *res_name, const char *atom_name)
{
    assert(classes); assert(res_name); assert(atom_name);
    int res, atom, status;
    status = freesasa_classify_find_atom(classes,res_name,atom_name,&res,&atom);
    if (status == FREESASA_SUCCESS)
        return classes->atom_class[res][atom];
    return freesasa_fail("%s: couldn't find classification of atom '%s %s'.",
                         __func__, res_name, atom_name);
}

const char* freesasa_classify_user_class2str(const freesasa_classify *classes, 
                                               int class)
{
    assert(classes);
    assert(class >= 0 && class < classes->n_classes);
    return classes->classes[class];
}
