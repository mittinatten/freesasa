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
		      const char chain_label,
		      const double x, 
		      const double y, 
		      const double z);

/** get array of coordinates */
const vector3* protein_xyz(const protein *p);

/** get array of radii */
const double* protein_r(const protein *p);

/** get number of atoms */
size_t protein_n(const protein *p);





#endif
