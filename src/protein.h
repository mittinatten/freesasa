#ifndef PROTEIN_H
#define PROTEIN_H

#include "vector3.h"

typedef struct protein protein;
typedef struct atom atom; 

/** init empty protein */
protein* protein_init();

/** free allocated memory */
void protein_free(protein*);

/** init protein with coordinates from pdb-file, 
    automatically skips hydrogens */
protein* protein_init_from_pdb(FILE *pdb_file);

/** string residue numbers allow for nonstandard formats, will include
    hydrogens if added (i.e. up to caller to make sure these are
    excluded if necessesary) */
void protein_add_atom(protein *p, 
		      const char* atom_name,
		      const char* residue_name, 
		      const char* residue_number,
		      char chain_label,
		      double x, double y, double z);

/** not tested yet */
void protein_update_atom(protein *p, int number, vector3 *coord);

/** not tested yet */
void protein_update_atom_xyz(protein *p, int number, 
			 double x, double y, double z);
/** not tested yet */
void protein_update_atoms(protein *p, vector3**);

/** get array of coordinates */
const vector3* protein_xyz(const protein *p);

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void protein_r(const protein *p, 
	       double *r,
	       double (*atom2radius)(const char *res_name, 
				     const char *atom_name));

/** get array of radii using default OONS radii, the array r is
    assumed to of correct size already */
void protein_r_def(const protein *p, double *r);

/** get number of atoms */
int protein_n(const protein *p);

/** get array of atoms */
const atom* protein_atoms(const protein *p);

/** get name of atom i */
const char* protein_atom_name(const protein *p, int i);

/** get name of residue atom i belongs to */
const char* protein_atom_res_name(const protein *p, int i);

/** get number of residue atom i belongs to */
const char* protein_atom_res_number(const protein *p, int i);

/** get chain atom i belongs to */
char protein_atom_chain(const protein *p, int i);

#endif
