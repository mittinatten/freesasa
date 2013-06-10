#ifndef ATOM_H
#define ATOM_H

#define ATOM_NAME_STRL 4
#define ATOM_RES_NAME_STRL 3
#define ATOM_RES_NUMBER_STRL 4

#include "vector3.h"

typedef enum {
    hydrogen, aliphatic_C, aromatic_C,
    carbo_C, amide_N, carbo_O, 
    hydroxyl_O, sulfur, selenium,
    unknown_polar, atom_type_unknown
} atom_type;

typedef enum {
    apolar, polar, charged, atom_class_unknown
} atom_class;

typedef struct {
    atom_type at;
    atom_class ac;
    char res_name[ATOM_RES_NAME_STRL];
    char res_number[ATOM_RES_NUMBER_STRL];
    char atom_name[ATOM_NAME_STRL];
    char chain_label;
    vector3 *xyz; //coordinates should be stored elsewhere (for speed)
} atom;

double atom_type2radius(atom_type);

atom_class atom_type2class(atom_type); 

/* Give the padded strings as arguments, i.e. positions 13-16 and
   18-20 of an PDB ATOM entry as atom_name and res_name */
atom_type atom_name2type(const char* res_name, const char* atom_name);

#endif
