/*
  Copyright Simon Mitternacht 2013-2015.

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

/**
   @file 
   @author Simon Mitternacht

   This header contains some functions to deal with
   ::freesasa_structure objects that for various reasons are not part
   of the public API.
 */

#include <stdio.h>
#include "freesasa.h"
#include "coord.h"

//! The maximum string length of returned by freesasa_structure_descriptor() 
#define FREESASA_STRUCTURE_DESCRIPTOR_STRL 20

/**
    Get coordinates.
    
    @param s Self.
    @return The coordinates of the structure as a ::freesasa_coord struct.
 */
const freesasa_coord* freesasa_structure_xyz(const freesasa_structure *s);

/**
    Get number of residues.

    Calculated crudely by determining the Number of unique
    combinations of residue name and chain label contained in the
    structure. If residues are mingled i.e. atoms of the same residue
    are present more than once at different places in the file this
    might be off.

    @param s Self.
    @return Number of residues.

    @ingroup StructureAPI
 */
int freesasa_structure_n_residues(const freesasa_structure *s);

/**
    Get a string describing an atom. 
    Format: "A    1 ALA  CA " 
    (chain label, residue number, residue type, atom name)

    @param s Self.
    @param i Atom index
    @return Descriptor string. 
 */
const char* freesasa_structure_atom_descriptor(const freesasa_structure *s, int i);

/**
    Get indices of first and last atoms of a residue
 
    @param s Self.
    @param r_i Residue index.
    @param first First atom of residue `r_i` will be stored here.
    @param last Last atom of residue `r_i` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if index `r_i` is invalid.
 */
int freesasa_structure_residue_atoms(const freesasa_structure *s, int r_i, 
                                     int *first, int *last);

/**
    Get a string describin a residue.
    Format: "A    1 ALA" (chain lable, resiude number, atom name)
    
    @param s Self.
    @param r_i atom index
    @return Descriptor string
 */
const char* freesasa_structure_residue_descriptor(const freesasa_structure *s, int r_i);

/**
    Writes PDB file, but with B-factors replaced by new values. Can be
    used to visualize SASA.

    @param s Self.
    @param output Output file.
    @param values Array of values to use as "B-factors" in the output.
    @param radii Array of atomic radii to put in the occupancy field in the output.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if there are problems with
    output-file.
 */
int freesasa_structure_write_pdb_bfactors(const freesasa_structure *s,
                                          FILE *output,
                                          const double *values,
                                          const double *radii);

#endif
