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
    double *atom_type_radius; // radius of each atom class
    double **atom_radius; // radius of each atom in each residue
};

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
    for (int i = 0; i < array_size; ++i) {
        assert(array[i]);
        if (strcmp(array[i],key) == 0) return i;
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

double freesasa_classify_user_radius(const freesasa_classify *classes, 
                                       const char *res_name, const char *atom_name)
{
    assert(classes); assert(res_name); assert(atom_name);
    int res = find_string(classes->residues,res_name,classes->n_residues);
    if (res < 0) {
        freesasa_warn("Unknown residue '%s', cannot calculate atom radius.", res_name);
        return -1.0;
    }
    int atom = find_string(classes->atoms[res],atom_name,classes->n_atoms[res]);
    if (atom < 0) {
        freesasa_warn("Unknown combination of residue and atom '%s %s', "
                      "cannot calculate atom radius.", res_name, atom_name);
        return -1.0;
    }
    return classes->atom_radius[res][atom];
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
    int res = find_string(classes->residues,res_name,classes->n_residues);
    if (res < 0) {
        freesasa_warn("Unknown residue '%s', cannot calculate atom radius.", res_name);
        return -1;
    }
    int atom = find_string(classes->atoms[res],atom_name,classes->n_atoms[res]);
    if (atom < 0) {
        freesasa_warn("Unknown combination of residue and atom '%s %s', "
                      "cannot calculate atom radius.", res_name, atom_name);
        return -1;
    }
    return classes->atom_class[res][atom];
}

const char* freesasa_classify_user_class2str(const freesasa_classify *classes, 
                                               int class)
{
    assert(classes);
    assert(class >= 0 && class < classes->n_classes);
    return classes->classes[class];
}
