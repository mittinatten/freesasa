#ifndef OONS_H
#define OONS_H

#include "atomclassifier.h"

double oons_radius_pdbline(const char *pdb_line);

double oons_radius(const char *res_name, const char *atom_name);

/** polar/apolar/unknown */
atomclassifier oons_classes();

/** carbo_O/aliphatic_C/aromatic_C/etc */
atomclassifier oons_types();


#endif
