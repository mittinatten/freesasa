#ifndef PROTEIN_H
#define PROTEIN_H

#include "vector3.h"

typedef struct protein protein;
typedef struct atom atom; 

/** init empty protein */
protein* protein_init();

/** free allocated memory */
void protein_free(protein*);

/** init protein with coordinates from pdb-file */
protein* protein_init_from_pdb(FILE *pdb_file);

/** string residue numbers allow for nonstandard formats */
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

/** get array of radii */
const double* protein_r(const protein *p);

/** get number of atoms */
size_t protein_n(const protein *p);

/** get array of atoms */
const atom* protein_atoms(const protein *p);



#endif
