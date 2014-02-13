/*
  Copyright Simon Mitternacht 2013-2014.

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

#include <stdio.h>
#include "coord.h"

typedef struct sasalib_structure_ sasalib_structure_t;

/** init empty protein */
sasalib_structure_t* sasalib_structure_init();

/** free allocated memory */
void sasalib_structure_free(sasalib_structure_t*);

/** Init protein with coordinates from pdb-file.  Automatically skips
    hydrogens. If an atom has alternative coordinates, only the first
    alternative is used. If a file has more than one MODEL (as in NMR
    structures) only the first model is used. If non-default behavior
    is wanted the pdb-file needs to be modified before calling this
    function, or atoms can be added manually using
    sasalib_structure_add_atom(). Returns NULL and prints error if
    input is invalid. */
sasalib_structure_t* sasalib_structure_init_from_pdb(FILE *pdb_file);

/** storing residue numbers as strings allows for nonstandard formats,
    will include hydrogens if added (i.e. up to caller to make sure
    these are excluded if necessesary) */
void sasalib_structure_add_atom(sasalib_structure_t *p, 
				const char* atom_name,
				const char* residue_name, 
				const char* residue_number,
				char chain_label,
				double x, double y, double z);

/** get array of coordinates */
const sasalib_coord_t* sasalib_structure_xyz(const sasalib_structure_t *p);

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void sasalib_structure_r(double *r,
			 const sasalib_structure_t *p, 
			 double (*atom2radius)(const char *res_name, 
					       const char *atom_name));

/** get array of radii using classify_radius(), the array r is assumed
    to of correct size already */
void sasalib_structure_r_def(double *r, const sasalib_structure_t *p);

/** get number of atoms */
int sasalib_structure_n(const sasalib_structure_t *p);

/** get name of atom i */
const char* sasalib_structure_atom_name(const sasalib_structure_t *p, 
					int i);

/** get name of residue atom i belongs to */
const char* sasalib_structure_atom_res_name(const sasalib_structure_t *p, 
					    int i);

/** get number of residue atom i belongs to */
const char* sasalib_structure_atom_res_number(const sasalib_structure_t *p, 
					      int i);

/** get chain atom i belongs to */
char sasalib_structure_atom_chain(const sasalib_structure_t *p, int i);

/** Writes PDB file based on structure, but with B-factors replaced by
    new values. Returns SASALIB_FAIL if there are problems with
    output-file. */
int sasalib_structure_write_pdb_bfactors(FILE *output, 
					 const sasalib_structure_t *p,
					 const double *values);

#endif
