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

   This header defines the type ::freesasa_structure, which
   represents a protein structure, and functions to deal with it.
 */

#include <stdio.h>
#include "coord.h"

//! The maximum string length of returned by freesasa_structure_descriptor() 
#define FREESASA_STRUCTURE_DESCRIPTOR_STRL 20

/**
    Struct for structure object.

    The struct includes coordinates, and atom names, etc. If it was
    initiated from a PDB file enough info will be stored so that
    a new PDB-file can be printed.
*/
typedef struct freesasa_structure freesasa_structure;

/**
    Allocate and initialize empty structure.

    @return The generated struct.
 */
freesasa_structure* freesasa_structure_new(void);

/**
    Free structure.
    
    @param s Self.
 */
void freesasa_structure_free(freesasa_structure *s);

/**
    Init protein with coordinates from pdb-file.  

    Reads in a PDB-file and generates a structure
    object. Automatically skips hydrogens. If an atom has alternative
    coordinates, only the first alternative is used. If a file has
    more than one `MODEL` (as in NMR structures) only the first model
    is used. User specifies if `HETATM` entries should be included. If
    non-default behavior is wanted, the PDB-file needs to be modified
    before calling this function, or atoms can be added manually using
    freesasa_structure_add_atom().

    @param pdb_file Input PDB-file.
    @param include_hetatm The value 0 means only read `ATOM` entries, 1 
    means also include `HETATM` entries.
    @return The generated struct. Returns `NULL` and prints error if
    input is invalid.
*/
freesasa_structure* freesasa_structure_from_pdb(FILE *pdb_file,
                                                int include_hetatm);

/**
    Add individual atom to structure.
    
    A structures can be built by adding atoms one by one. Storing
    residue numbers as strings allows for non-numeric labels. Will
    include hydrogens if added (i.e. up to caller to make sure these
    are excluded if necessesary).

    @param s Self.
    @param atom_name String of the format `" CA "`, `" OXT"`, etc.
    @param residue_name String of the format `"ALA"`, `"PHE"`, etc.
    @param residue_number String of the format `"   1"`, `" 123"`, etc.
    @param chain_label Any character to label chain, typically `'A'`, `'B'`, etc.
    @param x x-coordinate of atom.
    @param y y-coordinate of atom.
    @param z z-coordinate of atom.
    @return ::FREESASA_SUCCESS if input valid. ::FREESASA_FAIL if any of
    the strings are malformatted. ::FREESASA_WARN if the atom type is
    unknown. */
int freesasa_structure_add_atom(freesasa_structure *s,
                                const char* atom_name,
                                const char* residue_name,
                                const char* residue_number,
                                char chain_label,
                                double x, double y, double z);

/**
    Get coordinates.
    
    @param s Self.
    @return The coordinates of the structure as a ::freesasa_coord struct.
 */
const freesasa_coord* freesasa_structure_xyz(const freesasa_structure *s);

/**
    Get array of radii using custom conversion function.     
    
    @param r The array were the results will be stored. 
    @param s Self.
    @param atom2radius Function pointer to function that converts 
      an atom label into a radius.
    @see freesasa_structure_r_def()
*/
void freesasa_structure_r(double *r,
                          const freesasa_structure *s,
                          double (*atom2radius)(const char *res_name,
                                                const char *atom_name));

/**
    Get array of default atomic radii. 

    Calls freesasa_structure_r() using freesasa_classify_radius() for
    radius calculations..
    
    @param r The array were the results will be stored.
    @param s Self.
    @see freesasa_classify_radius()
*/

void freesasa_structure_r_def(double *r, const freesasa_structure *s);

/**
    Get number of atoms.
    
    @param s Self.
    @return Number of atoms.
*/
int freesasa_structure_n(const freesasa_structure *s);

/**
    Get number of residues.

    Number of unique combinations of residue name and chain label
    contained in the structure.

    @param s Self.
    @return Number of residues.
 */
int freesasa_structure_n_residues(const freesasa_structure *s);
/**
    Get atom name
    
    @param s Self.
    @param i Atom index.
    @return Atom name in the form `" CA "`, `" OXT"`, etc.
 */
const char* freesasa_structure_atom_name(const freesasa_structure *s,
                                         int i);

/**
    Get residue name.
    
    @param s Self.
    @param i Atom index.
    @return Residue name in the form `"ALA"`, `"PHE"`, etc.
*/
const char* freesasa_structure_atom_res_name(const freesasa_structure *s,
                                             int i);

/**
    Get residue number.

    @param s Self.
    @param i Atom index.
    @return Residue name in the form `"   1"`, `" 123"`, etc.
*/
const char* freesasa_structure_atom_res_number(const freesasa_structure *s,
                                               int i);

/**
    Get chain label.
   
    @param s Self.
    @param i Atom index.
    @return Chain label (`'A'`, `'B'`, etc.)
 */
char freesasa_structure_atom_chain(const freesasa_structure *s, int i);

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
