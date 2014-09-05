/*
  Copyright Simon Mitternacht 2013-2014.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef FREESASA_STRUCTURE_H
#define FREESASA_STRUCTURE_H

#include <stdio.h>
#include "coord.h"

typedef struct freesasa_structure_ freesasa_structure_t;

/** init empty protein */
freesasa_structure_t* freesasa_structure_init();

/** free allocated memory */
void freesasa_structure_free(freesasa_structure_t*);

/** Init protein with coordinates from pdb-file.  Automatically skips
    hydrogens. If an atom has alternative coordinates, only the first
    alternative is used. If a file has more than one MODEL (as in NMR
    structures) only the first model is used. If non-default behavior
    is wanted the pdb-file needs to be modified before calling this
    function, or atoms can be added manually using
    freesasa_structure_add_atom(). Returns NULL and prints error if
    input is invalid. */
freesasa_structure_t* freesasa_structure_init_from_pdb(FILE *pdb_file);

/** storing residue numbers as strings allows for nonstandard formats,
    will include hydrogens if added (i.e. up to caller to make sure
    these are excluded if necessesary) */
void freesasa_structure_add_atom(freesasa_structure_t *p,
                                 const char* atom_name,
                                 const char* residue_name,
                                 const char* residue_number,
                                 char chain_label,
                                 double x, double y, double z);

/** get array of coordinates */
const freesasa_coord_t* freesasa_structure_xyz(const freesasa_structure_t *p);

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void freesasa_structure_r(double *r,
                          const freesasa_structure_t *p,
                          double (*atom2radius)(const char *res_name,
                                                const char *atom_name));

/** get array of radii using classify_radius(), the array r is assumed
    to of correct size already */
void freesasa_structure_r_def(double *r, const freesasa_structure_t *p);

/** get number of atoms */
int freesasa_structure_n(const freesasa_structure_t *p);

/** get name of atom i */
const char* freesasa_structure_atom_name(const freesasa_structure_t *p,
                                         int i);

/** get name of residue atom i belongs to */
const char* freesasa_structure_atom_res_name(const freesasa_structure_t *p,
                                             int i);

/** get number of residue atom i belongs to */
const char* freesasa_structure_atom_res_number(const freesasa_structure_t *p,
                                               int i);

/** get chain atom i belongs to */
char freesasa_structure_atom_chain(const freesasa_structure_t *p, int i);

/** Writes PDB file based on structure, but with B-factors replaced by
    new values. Returns FREESASA_FAIL if there are problems with
    output-file. */
int freesasa_structure_write_pdb_bfactors(FILE *output,
                                          const freesasa_structure_t *p,
                                          const double *values);

#endif
