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
#include "protein.h"
#include "pdbutil.h"
#include "oons.h"


struct atom {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char chain_label;
};

struct protein {
    atom *a;
    vector3 *xyz;
    int number_atoms;
    int number_residues;
    int number_chains;
};

protein* protein_init()
{
    protein *p = (protein*) malloc(sizeof(protein));
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = (atom*) malloc(sizeof(atom));
    p->xyz = (vector3*) malloc(sizeof(vector3));
    return p;
}

void protein_free(protein *p)
{
    free(p->a);
    free(p->xyz);
    free(p);
}

/* returns alt_label if there is any */
char protein_get_pdb_atom(const char *line, atom *a, vector3 *v)
{
    a->chain_label = pdbutil_get_chain_label(line);
    pdbutil_get_coord(line, v);
    pdbutil_get_atom_name(line, a->atom_name);
    pdbutil_get_res_name(line, a->res_name);
    pdbutil_get_res_number(line, a->res_number);
    return pdbutil_get_alt_coord_label(line);
}

protein* protein_init_from_pdb(FILE *pdb_file)
{
    protein *p = protein_init();
    /* two possible implementations, either read file in two passes,
       first to determine number of atoms, second to read them in,
       keeping the number of mallocs/reallocs to a minimum.  Or read
       file once using protein_add_atom() for realloc.  Will begin
       with second alternative since that has to be implemented
       anyway. */
    size_t len = 80;
    char *line = (char*) malloc(sizeof(char)*(len+1));
    char the_alt = ' ';
    while (getline(&line, &len, pdb_file) != -1) {
        if (strncmp("ATOM",line,4)==0) {
	    if (pdbutil_ishydrogen(line)) continue;
            vector3 v;
            atom a;
	    char alt = protein_get_pdb_atom(line,&a,&v);
	    if ((alt != ' ' && the_alt == ' ') || (alt == ' ')) { 
		the_alt = alt;
	    } else if (alt != ' ' && alt != the_alt) {
		continue;
	    } 
	    
	    protein_add_atom(p,a.atom_name,a.res_name,a.res_number,
			     a.chain_label, v.x, v.y, v.z);
	}
	if (strncmp("ENDMDL",line,4)==0) break;
    }
    free(line);
    return p;
}

void protein_alloc_one(protein *p)
{
    int na = ++p->number_atoms;
    p->a = (atom*) realloc(p->a,sizeof(atom)*na);
    p->xyz = (vector3*) realloc(p->xyz,sizeof(vector3)*na);
}

void protein_add_atom(protein *p,
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
    protein_alloc_one(p);

    int na = p->number_atoms;

    vector3_set(&p->xyz[na-1],x,y,z);

    atom *a = &p->a[na-1];
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
void protein_r(double *r,
	       const protein *p, 
	       double (*atom2radius)(const char *res_name, 
				     const char *atom_name))
{
    for (int i = 0; i < p->number_atoms; ++i) {
	r[i] = atom2radius(p->a[i].res_name, p->a[i].atom_name);
    }
}

/** get array of radii using default OONS radii, the array r is
    assumed to of correct size already */
void protein_r_def(double *r, const protein *p)
{
    protein_r(r,p,oons_radius);
}

void protein_update_atom_xyz(protein* p, int number,
			     double x, double y, double z)
{
    assert(number < p->number_atoms);
    vector3_set(&p->xyz[number], x, y, z);
}

void protein_update_atom(protein* p, int number, vector3 *v)
{
    assert(number < p->number_atoms);
    vector3_copy(&p->xyz[number], v);
}

void protein_update_atoms(protein* p, vector3 **v) 
{
    for (int i = 0; i < p->number_atoms; ++i) {
	vector3_copy(&p->xyz[i],v[i]);
    }
}

const vector3* protein_xyz(const protein *p)
{
    return p->xyz;
}

int protein_n(const protein *p)
{
    return p->number_atoms;
}

const char* protein_atom_name(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].atom_name;
}

const char* protein_atom_res_name(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_name;
}

const char* protein_atom_res_number(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_number;
}

char protein_atom_chain(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].chain_label;
}
