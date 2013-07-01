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

double oons_radius_pdbline(const char *pdb_line);

double oons_radius(const char *res_name, const char *atom_name);

/** polar/apolar/unknown */
atomclassifier oons_classes();

/** carbo_O/aliphatic_C/aromatic_C/etc */
atomclassifier oons_types();


#endif
