/*
  Copyright Simon Mitternacht 2013.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

// needed for getline()
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "structure.h"
#include "pdbutil.h"
#include "classify.h"
#include "coord.h"

typedef struct {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char chain_label;
} atom_t;

struct structure_ {
    atom_t *a;
    coord_t *xyz;
    int number_atoms;
    int number_residues;
    int number_chains;
};

structure_t* structure_init()
{
    structure_t *p = (structure_t*) malloc(sizeof(structure_t));
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = (atom_t*) malloc(sizeof(atom_t));
    p->xyz = coord_new();
    return p;
}

void structure_free(structure_t *p)
{
    free(p->a);
    coord_free(p->xyz);
    free(p);
}

/* returns alt_label if there is any */
char structure_get_pdb_atom(atom_t *a, double *xyz, const char *line)
{
    a->chain_label = pdbutil_get_chain_label(line);
    pdbutil_get_coord(xyz, line);
    pdbutil_get_atom_name(a->atom_name, line);
    pdbutil_get_res_name(a->res_name, line);
    pdbutil_get_res_number(a->res_number, line);
    return pdbutil_get_alt_coord_label(line);
}

structure_t* structure_init_from_pdb(FILE *pdb_file)
{
    structure_t *p = structure_init();
    /* two possible implementations, either read file in two passes,
       first to determine number of atoms, second to read them in,
       keeping the number of mallocs/reallocs to a minimum.  Or read
       file once using structure_add_atom() for realloc.  Will begin
       with second alternative since that has to be implemented
       anyway. */
    size_t len = 80;
    char *line = (char*) malloc(sizeof(char)*(len+1));
    char the_alt = ' ';
    while (getline(&line, &len, pdb_file) != -1) {
        if (strncmp("ATOM",line,4)==0) {
            if (pdbutil_ishydrogen(line)) continue;
            double v[3];
            atom_t a;
            char alt = structure_get_pdb_atom(&a,v,line);
            if ((alt != ' ' && the_alt == ' ') || (alt == ' ')) { 
                the_alt = alt;
            } else if (alt != ' ' && alt != the_alt) {
                continue;
            } 
            
            structure_add_atom(p,a.atom_name,a.res_name,a.res_number,
                               a.chain_label, v[0], v[1], v[2]);
        }
        if (strncmp("ENDMDL",line,4)==0) break;
    }
    free(line);
    return p;
}

static void structure_alloc_one(structure_t *p)
{
    int na = ++p->number_atoms;
    p->a = (atom_t*) realloc(p->a,sizeof(atom_t)*na);
}

void structure_add_atom(structure_t *p,
                        const char *atom_name,
                        const char *residue_name,
                        const char *residue_number,
                        char chain_label,
                        double x, double y, double z)
{
    // check input for consistency
    assert(strlen(atom_name) == PDB_ATOM_NAME_STRL);
    assert(strlen(residue_name) == PDB_ATOM_RES_NAME_STRL);
    assert(strlen(residue_number) == PDB_ATOM_RES_NUMBER_STRL);

    // allocate memory, increase number of atoms counter
    structure_alloc_one(p);
    int na = p->number_atoms;
    coord_append_xyz(p->xyz,&x,&y,&z,1);
    atom_t *a = &p->a[na-1];
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    a->chain_label = chain_label;

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

    // what else?
    // deal with alternate location indicators...
}

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void structure_r(double *r,
		 const structure_t *p, 
		 double (*atom2radius)(const char *res_name, 
				       const char *atom_name))
{
    for (int i = 0; i < p->number_atoms; ++i) {
	r[i] = atom2radius(p->a[i].res_name, p->a[i].atom_name);
    }
}

/** get array of radii using default OONS radii, the array r is
    assumed to of correct size already */
void structure_r_def(double *r, const structure_t *p)
{
    structure_r(r,p,classify_radius);
}

const coord_t* structure_xyz(const structure_t *p)
{
    return p->xyz;
}

int structure_n(const structure_t *p)
{
    return p->number_atoms;
}

const char* structure_atom_name(const structure_t *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].atom_name;
}

const char* structure_atom_res_name(const structure_t *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_name;
}

const char* structure_atom_res_number(const structure_t *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_number;
}

char structure_atom_chain(const structure_t *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].chain_label;
}
