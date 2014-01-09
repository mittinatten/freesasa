/*
  Copyright Simon Mitternacht 2013-2014.

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
#include <errno.h>
#include "structure.h"
#include "pdb.h"
#include "classify.h"
#include "coord.h"

extern const char *sasalib_name;

typedef struct {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char line[PDB_LINE_STRL+1];
    char chain_label;
} atom_t;

struct sasalib_structure_ {
    atom_t *a;
    sasalib_coord_t *xyz;
    int number_atoms;
    int number_residues;
    int number_chains;
};

sasalib_structure_t* sasalib_structure_init()
{
    sasalib_structure_t *p 
	= (sasalib_structure_t*) malloc(sizeof(sasalib_structure_t));
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = (atom_t*) malloc(sizeof(atom_t));
    p->xyz = sasalib_coord_new();
    return p;
}

void sasalib_structure_free(sasalib_structure_t *p)
{
    free(p->a);
    sasalib_coord_free(p->xyz);
    free(p);
}

/* returns alt_label if there is any */
char sasalib_structure_get_pdb_atom(atom_t *a, double *xyz, const char *line)
{
    a->chain_label = sasalib_pdb_get_chain_label(line);
    sasalib_pdb_get_coord(xyz, line);
    sasalib_pdb_get_atom_name(a->atom_name, line);
    sasalib_pdb_get_res_name(a->res_name, line);
    sasalib_pdb_get_res_number(a->res_number, line);
    return sasalib_pdb_get_alt_coord_label(line);
}

sasalib_structure_t* sasalib_structure_init_from_pdb(FILE *pdb_file)
{
    sasalib_structure_t *p = sasalib_structure_init();

    size_t len = PDB_LINE_STRL;
    char *line = (char*) malloc(sizeof(char)*(len+1));
    char the_alt = ' ';
    while (getline(&line, &len, pdb_file) != -1) {
        if (strncmp("ATOM",line,4)==0) {
            if (sasalib_pdb_ishydrogen(line)) continue;
            double v[3];
            atom_t a;
            char alt = sasalib_structure_get_pdb_atom(&a,v,line);
            if ((alt != ' ' && the_alt == ' ') || (alt == ' ')) { 
                the_alt = alt;
            } else if (alt != ' ' && alt != the_alt) {
                continue;
            } 
            
            sasalib_structure_add_atom(p,a.atom_name,a.res_name,a.res_number,
				       a.chain_label, v[0], v[1], v[2]);
	    strncpy(p->a[p->number_atoms-1].line,line,PDB_LINE_STRL);
        }
        if (strncmp("ENDMDL",line,4)==0) {
	    if (p->number_atoms == 0) {
		fprintf(stderr, "%s: error: input had ENDMDL before "
			"first ATOM entry.\n", sasalib_name);
		free(line);
		sasalib_structure_free(p);
		return NULL;
	    }
	    break;
	}
    }
    free(line);
    if (p->number_atoms == 0) {
	fprintf(stderr,"%s: error: input had no ATOM entries.\n",
		sasalib_name);
	sasalib_structure_free(p);
	return NULL;
    }
    return p;
}

static void sasalib_structure_alloc_one(sasalib_structure_t *p)
{
    int na = ++p->number_atoms;
    p->a = (atom_t*) realloc(p->a,sizeof(atom_t)*na);
}

void sasalib_structure_add_atom(sasalib_structure_t *p,
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
    sasalib_structure_alloc_one(p);
    int na = p->number_atoms;
    sasalib_coord_append_xyz(p->xyz,&x,&y,&z,1);
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

void sasalib_structure_r(double *r,
			 const sasalib_structure_t *p, 
			 double (*atom2radius)(const char *res_name, 
					       const char *atom_name))
{
    for (int i = 0; i < p->number_atoms; ++i) {
	r[i] = atom2radius(p->a[i].res_name, p->a[i].atom_name);
    }
}

void sasalib_structure_r_def(double *r, const sasalib_structure_t *p)
{
    sasalib_structure_r(r,p,sasalib_classify_radius);
}

const sasalib_coord_t* sasalib_structure_xyz(const sasalib_structure_t *p)
{
    return p->xyz;
}

int sasalib_structure_n(const sasalib_structure_t *p)
{
    return p->number_atoms;
}

const char* sasalib_structure_atom_name(const sasalib_structure_t *p, 
					int i)
{
    assert (i < p->number_atoms);
    return p->a[i].atom_name;
}

const char* sasalib_structure_atom_res_name(const sasalib_structure_t *p, 
					    int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_name;
}

const char* sasalib_structure_atom_res_number(const sasalib_structure_t *p, 
					      int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_number;
}

char sasalib_structure_atom_chain(const sasalib_structure_t *p, 
				  int i)
{
    assert (i < p->number_atoms);
    return p->a[i].chain_label;
}

int sasalib_structure_write_pdb_bfactors(FILE *output, 
					 const sasalib_structure_t *p,
					 double *values)
{
    if (!output) {
	fprintf(stderr,
		"%s: error: NULL file pointer provided for PDB output.\n",
		sasalib_name);
	return SASALIB_FAIL;
    }
    // Write ATOM entries
    char buf[PDB_LINE_STRL+1];
    int n = sasalib_structure_n(p);
    for (int i = 0; i < n; ++i) {
	strncpy(buf,p->a[i].line,PDB_LINE_STRL);
	sprintf(&buf[60],"%6.2f",values[i]); 
	errno = 0;

	if (fprintf(output,"%s\n",buf)<0) {
	    fprintf(stderr,"%s: error: Problem writing new PDB-file. %s\n",
		    sasalib_name, strerror(errno));
	    return SASALIB_FAIL;
	}
    }
    // Write TER line
    errno = 0;
    char buf2[6];
    strncpy(buf2,&buf[6],5);
    buf2[5]='\0';
    if (fprintf(output,"TER   %5d     %4s %c%4s\n",
		atoi(buf2)+1, p->a[n-1].res_name,
		p->a[n-1].chain_label, p->a[n-1].res_number) < 0) {
	fprintf(stderr,"%s: error: Problem writing new PDB-file. %s\n",
		sasalib_name, strerror(errno));
	return SASALIB_FAIL;
    }

    return SASALIB_SUCCESS;
}
