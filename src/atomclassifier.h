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

#ifndef CLASSIFIER_H
#define CLASSIFIER_H

/** this struct can be used to pass information about how to classify atoms
    based on residue and the atom name (in PDB format).
 */

typedef struct {
    int nclasses;
    const char **class2str;
    const char *name;
    int (*classify)(const char *res_name, const char *atom_name);
} atomclassifier;

/** This is an example classifier for identifying residues (see
    atomclassifier.c for implementation). In addition to the 20
    standard amino acids it includes ASX, GLX, CSE/SEC and
    PYH/PYL. (CSE and PYH have not been tested, haven't found any
    example PDB files). */
atomclassifier atomclassifier_residue();

#endif
