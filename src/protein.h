#ifndef PROTEIN_H
#define PROTEIN_H

#include "atom.h"
#include "vector3.h"

typedef struct {
    atom *a;
    vector3 *xyz;
    double *r;
    size_t number_atoms;
    size_t number_residues;
    size_t number_chains;
} protein;


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

#endif
