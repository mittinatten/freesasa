#ifndef ATOM_H
#define ATOM_H

#include "vector3.h"

typedef enum {
	hydrogen, aliphatic_C, aromatic_C,
    carbo_C, amide_N, carbo_O, 
	hydroxyl_O, sulfur, selenium,
	unknown_polar, unknown
} atom_type;

typedef enum {
	apolar, polar, charged, unknown
} atom_class;

double atom_radius (atom_type);

atom_class atom_type2class (atom_type); 

/* Give the padded strings as arguments, i.e. positions 13-16 and
   18-20 of an PDB ATOM entry as atom_name and res_name */
atom_type atom_name2type (char* res_name, char* atom_name);

#endif
