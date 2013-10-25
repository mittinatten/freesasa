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

#ifndef SASALIB_STRUCTURE_H
#define SASALIB_STRUCTURE_H

#include "vector3.h"

typedef struct structure_ structure_t;

/** init empty protein */
structure_t* structure_init();

/** free allocated memory */
void structure_free(structure_t*);

/** Init protein with coordinates from pdb-file.  Automatically skips
    hydrogens. If an atom has alternative coordinates, only the first
    alternative is used. If a file has more than one MODEL (as in NMR
    structures) only the first model is used. If non-default behavior
    is wanted the pdb-file needs to be modified before calling this
    function, or atoms can be added manually using
    structure_add_atom(). */
structure_t* structure_init_from_pdb(FILE *pdb_file);

/** string residue numbers allow for nonstandard formats, will include
    hydrogens if added (i.e. up to caller to make sure these are
    excluded if necessesary) */
void structure_add_atom(structure_t *p, 
		      const char* atom_name,
		      const char* residue_name, 
		      const char* residue_number,
		      char chain_label,
		      double x, double y, double z);

/** not tested yet */
void structure_update_atom(structure_t *p, int number, vector3 *coord);

/** not tested yet */
void structure_update_atom_xyz(structure_t *p, int number, 
			 double x, double y, double z);
/** not tested yet */
void structure_update_atoms(structure_t *p, vector3**);

/** get array of coordinates */
const vector3* structure_xyz(const structure_t *p);

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void structure_r(double *r,
	       const structure_t *p, 
	       double (*atom2radius)(const char *res_name, 
				     const char *atom_name));

/** get array of radii using default OONS radii, the array r is
    assumed to of correct size already */
void structure_r_def(double *r, const structure_t *p);

/** get number of atoms */
int structure_n(const structure_t *p);

/** get name of atom i */
const char* structure_atom_name(const structure_t *p, int i);

/** get name of residue atom i belongs to */
const char* structure_atom_res_name(const structure_t *p, int i);

/** get number of residue atom i belongs to */
const char* structure_atom_res_number(const structure_t *p, int i);

/** get chain atom i belongs to */
char structure_atom_chain(const structure_t *p, int i);

#endif
