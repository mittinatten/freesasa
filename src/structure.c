/*
  Copyright Simon Mitternacht 2013-2014.

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

// needed for getline()
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "structure.h"
#include "pdb.h"
#include "classify.h"
#include "coord.h"

extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

typedef struct {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char line[PDB_LINE_STRL+1];
    char chain_label;
} atom_t;

struct freesasa_structure_ {
    atom_t *a;
    freesasa_coord_t *xyz;
    int number_atoms;
    int number_residues;
    int number_chains;
};

freesasa_structure_t* freesasa_structure_init()
{
    freesasa_structure_t *p
        = (freesasa_structure_t*) malloc(sizeof(freesasa_structure_t));
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = (atom_t*) malloc(sizeof(atom_t));
    p->xyz = freesasa_coord_new();
    return p;
}

void freesasa_structure_free(freesasa_structure_t *p)
{
    if (p == NULL) return;
    if (p->a) free(p->a);
    if (p->xyz) freesasa_coord_free(p->xyz);
    free(p);
}

/* returns alt_label if there is any */
char freesasa_structure_get_pdb_atom(atom_t *a, double *xyz, const char *line)
{
    a->chain_label = freesasa_pdb_get_chain_label(line);
    freesasa_pdb_get_coord(xyz, line);
    freesasa_pdb_get_atom_name(a->atom_name, line);
    freesasa_pdb_get_res_name(a->res_name, line);
    freesasa_pdb_get_res_number(a->res_number, line);
    return freesasa_pdb_get_alt_coord_label(line);
}

freesasa_structure_t* freesasa_structure_init_from_pdb(FILE *pdb_file)
{
    freesasa_structure_t *p = freesasa_structure_init();

    size_t len = PDB_LINE_STRL;
    char *line = (char*) malloc(sizeof(char)*(len+1));
    char the_alt = ' ';
    while (getline(&line, &len, pdb_file) != -1) {
        if (strncmp("ATOM",line,4)==0) {
            if (freesasa_pdb_ishydrogen(line)) continue;
            double v[3];
            atom_t a;
            char alt = freesasa_structure_get_pdb_atom(&a,v,line);
            if ((alt != ' ' && the_alt == ' ') || (alt == ' ')) {
                the_alt = alt;
            } else if (alt != ' ' && alt != the_alt) {
                continue;
            }

            freesasa_structure_add_atom(p,a.atom_name,a.res_name,a.res_number,
                                        a.chain_label, v[0], v[1], v[2]);
            strncpy(p->a[p->number_atoms-1].line,line,PDB_LINE_STRL);
        }
        if (strncmp("ENDMDL",line,4)==0) {
            if (p->number_atoms == 0) {
                freesasa_fail("input had ENDMDL before first ATOM entry.");
                free(line);
                freesasa_structure_free(p);
                return NULL;
            }
            break;
        }
    }
    free(line);
    if (p->number_atoms == 0) {
        freesasa_fail("input had no ATOM entries.");
        freesasa_structure_free(p);
        return NULL;
    }
    return p;
}

static void freesasa_structure_alloc_one(freesasa_structure_t *p)
{
    int na = ++p->number_atoms;
    p->a = (atom_t*) realloc(p->a,sizeof(atom_t)*na);
}

int freesasa_structure_add_atom(freesasa_structure_t *p,
                                const char *atom_name,
                                const char *residue_name,
                                const char *residue_number,
                                char chain_label,
                                double x, double y, double z)
{
    // check input for consistency
    int validity = freesasa_classify_validate_atom(residue_name,atom_name);
    if (validity != FREESASA_SUCCESS) {
        return freesasa_warn("Skipping atom '%s' in residue '%s'",
                             atom_name,residue_name);
    }

    // allocate memory, increase number of atoms counter
    freesasa_structure_alloc_one(p);
    int na = p->number_atoms;
    freesasa_coord_append_xyz(p->xyz,&x,&y,&z,1);
    atom_t *a = &p->a[na-1];
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    a->chain_label = chain_label;
    a->line[0] = '\0';

    /* here we assume atoms are ordered sequentially, i.e. chains are
       not mixed in input: if two sequential atoms have different
       residue numbers or chain labels a new residue or new chain is
       assumed to begin */
    if (p->number_chains == 0)
        ++p->number_chains;
    if (na > 1 && chain_label != p->a[na-2].chain_label)
        ++p->number_chains;

    if (p->number_residues == 0)
        ++p->number_residues;
    if (na > 1 && strcmp(residue_number,p->a[na-2].res_number))
        ++p->number_residues;
}

void freesasa_structure_r(double *r,
                          const freesasa_structure_t *p,
                          double (*atom2radius)(const char *res_name,
                                                const char *atom_name))
{
    for (int i = 0; i < p->number_atoms; ++i) {
        r[i] = atom2radius(p->a[i].res_name, p->a[i].atom_name);
    }
}

void freesasa_structure_r_def(double *r, const freesasa_structure_t *p)
{
    freesasa_structure_r(r,p,freesasa_classify_radius);
}

const freesasa_coord_t* freesasa_structure_xyz(const freesasa_structure_t *p)
{
    return p->xyz;
}

int freesasa_structure_n(const freesasa_structure_t *p)
{
    return p->number_atoms;
}

const char* freesasa_structure_atom_name(const freesasa_structure_t *p,
                                         int i)
{
    assert (i < p->number_atoms);
    return p->a[i].atom_name;
}

const char* freesasa_structure_atom_res_name(const freesasa_structure_t *p,
                                             int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_name;
}

const char* freesasa_structure_atom_res_number(const freesasa_structure_t *p,
                                               int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_number;
}

char freesasa_structure_atom_chain(const freesasa_structure_t *p,
                                   int i)
{
    assert (i < p->number_atoms);
    return p->a[i].chain_label;
}

int freesasa_structure_write_pdb_bfactors(FILE *output,
                                          const freesasa_structure_t *p,
                                          const double *values)
{
    if (!output)
        return freesasa_fail("NULL file pointer provided for PDB output.");
    if (!values)
        return freesasa_fail("NULL file pointer provided for b-factors.");
    // Write ATOM entries
    char buf[PDB_LINE_STRL+1];
    int n = freesasa_structure_n(p);
    for (int i = 0; i < n; ++i) {
        if (p->a[i].line[0] == '\0') {
            return freesasa_fail("Attempting to write B-factors for "
                                 "structure not initialized from PDB-file. "
                                 "Aborting.");
        }
        strncpy(buf,p->a[i].line,PDB_LINE_STRL);
        sprintf(&buf[60],"%6.2f",values[i]);
        errno = 0;
        if (fprintf(output,"%s\n",buf)<0)
            return freesasa_fail("Problem writing new PDB-file. %s",
                                 strerror(errno));
    }
    // Write TER line
    errno = 0;
    char buf2[6];
    strncpy(buf2,&buf[6],5);
    buf2[5]='\0';
    if (fprintf(output,"TER   %5d     %4s %c%4s\n",
                atoi(buf2)+1, p->a[n-1].res_name,
                p->a[n-1].chain_label, p->a[n-1].res_number) < 0) {
        freesasa_fail("Problem writing new PDB-file. %s",
                      strerror(errno));
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;
}
