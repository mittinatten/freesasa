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

#ifndef OONS_H
#define OONS_H

#include "atomclassifier.h"

/** the return values of the atomclassifier-function in
    oons_classes() */
typedef enum {
    apolar=0, polar, oons_class_unknown
} oons_class;

/** the return values of the atomclassifier-function in
    oons_types() */
typedef enum {
    hydrogen=0, aliphatic_C, aromatic_C,
    carbo_C, amide_N, carbo_O, 
    hydroxyl_O, sulfur, selenium,
    unknown_polar, oons_type_unknown
} oons_type;

double oons_radius_pdbline(const char *pdb_line);

double oons_radius(const char *res_name, const char *atom_name);

/** polar/apolar/unknown */
atomclassifier oons_classes();

/** carbo_O/aliphatic_C/aromatic_C/etc */
atomclassifier oons_types();



#endif
