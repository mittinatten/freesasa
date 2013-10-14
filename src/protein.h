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

#ifndef PROTEIN_H
#define PROTEIN_H

#include "vector3.h"

typedef struct protein protein;
typedef struct atom atom; 

/** init empty protein */
protein* protein_init();

/** free allocated memory */
void protein_free(protein*);

/** Init protein with coordinates from pdb-file.  Automatically skips
    hydrogens. If an atom has alternative coordinates, only the first
    alternative is used. If a file has more than one MODEL (as in NMR
    structures) only the first model is used. If non-default behavior
    is wanted the pdb-file needs to be modified before calling this
    function, or atoms can be added manually using
    protein_add_atom(). */
protein* protein_init_from_pdb(FILE *pdb_file);

/** string residue numbers allow for nonstandard formats, will include
    hydrogens if added (i.e. up to caller to make sure these are
    excluded if necessesary) */
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

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void protein_r(double *r,
	       const protein *p, 
	       double (*atom2radius)(const char *res_name, 
				     const char *atom_name));

/** get array of radii using default OONS radii, the array r is
    assumed to of correct size already */
void protein_r_def(double *r, const protein *p);

/** get number of atoms */
int protein_n(const protein *p);

/** get name of atom i */
const char* protein_atom_name(const protein *p, int i);

/** get name of residue atom i belongs to */
const char* protein_atom_res_name(const protein *p, int i);

/** get number of residue atom i belongs to */
const char* protein_atom_res_number(const protein *p, int i);

/** get chain atom i belongs to */
char protein_atom_chain(const protein *p, int i);

#endif
