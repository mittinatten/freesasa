#ifndef PROTEIN_H
#define PROTEIN_H

#include "atom.h"
#include "vector3.h"

typedef struct {
	atom_type *at;
	atom_class *ac;
	char **res_name;
	char **atom_name;
	char **res_number;
	char *chain;
	vector3 *coord;
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
