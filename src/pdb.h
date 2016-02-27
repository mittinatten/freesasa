/*
  Copyright Simon Mitternacht 2013-2016.

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

#ifndef FREESASA_PDB_H
#define FREESASA_PDB_H

#include "freesasa.h"
#include "util.h"

/**
    @file
    @author Simon Mitternacht

    The following functions all extract info from the PDB lines `ATOM`
    and `HETATM`. Valid lines have to begin with either `ATOM` or
    `HETATM` and be sufficiently long to contain the value in
    question. 
*/


#define PDB_ATOM_NAME_STRL 4 //!< Length of strings with atom names, such as `" CA "`.
#define PDB_ATOM_RES_NAME_STRL 3 //!< Length of string with residue names, such as `"ALA"`.
#define PDB_ATOM_RES_NUMBER_STRL 4 //!< Length of string with residue number, such as `" 123"`.
#define PDB_ATOM_SYMBOL_STRL 2 //!< Length for string with element symbol, such "FE"
#define PDB_LINE_STRL 80 //!< Length of a line in PDB file.

/**
    Finds the location of all MODEL entries in the file pdb, returns
    the number of models found. 

    The array *ranges will be dynamically allocated to contain a
    file ranges for each model.

    @return Returns 0 if no MODEL lines were found (for example an
      X-ray structure with only one model) and sets *ranges to
      NULL. A return value of 0 doesn't have to mean the file is
      empty. ::FREESASA_FAIL if malloc-failure.
 */
int
freesasa_pdb_get_models(FILE* pdb,
                        struct file_range** ranges);

/**
    Finds the location of all chains within the file range 'model'.

    *ranges points to an array of dynamically allocated file
    ranges for each model.
    
 */
int
freesasa_pdb_get_chains(FILE *pdb,
                        struct file_range model,
                        struct file_range **ranges,
                        int options);
/**
    Get atom name from a PDB line.
    
    Extracts the whole atom name field from an `ATOM` or `HETATM` PDB
    line, including padding, i.e. a string of ::PDB_ATOM_NAME_STRL
    characters. For example `" CA "` for a regular C-alpha.  If the
    line is invalid, the name will be an empty string and the function
    returns ::FREESASA_FAIL.
    
    @param name The name is written to this string.
    @param line Line from a PDB file.
    @return ::FREESASA_SUCCESS if input is readable, else ::FREESASA_FAIL.
 */
int
freesasa_pdb_get_atom_name(char *name,
                           const char *line);

/**
    Get residue name from a PDB line.
   
    Extracts the whole residue name from an `ATOM` or `HETATM` PDB
    line, i.e. a string of ::PDB_ATOM_RES_NAME_STRL characters. For
    example `"ALA"` for an Alanine. If theline is invalid, the name
    will be an empty string and the function returns ::FREESASA_FAIL. 

    @param name The name is written to this string.
    @param line Line from a PDB file.
    @return ::FREESASA_SUCCESS if input is readable, else ::FREESASA_FAIL.
 */
int
freesasa_pdb_get_res_name(char *name,
                          const char *line);

/**
    Get atom coordinates from a PDB line.

    Extracts x-, y- and z-coordinates from an `ATOM` or `HETATM` PDB
    line. If the line is invalid, the function returns ::FREESASA_FAIL
    and coord will remain unchanged.

    @param coord The coordiantes are written to this array as x,y,z.
    @param line Line from a PDB file.
    @return ::FREESASA_SUCCESS if input is readable, else ::FREESASA_FAIL.
 */
int
freesasa_pdb_get_coord(double *coord,
                       const char *line);

/**
    Get residue number from a PDB line.
    
    Extracts residue number (ResSeq) as a string from an `ATOM` or
    `HETATM` PDB line as a string. String format is used because not
    all residue-numbers are numbers. The string should have length
    ::PDB_ATOM_RES_NUMBER_STRL. If the line is invalid, the number will be an
    empty string and the function returns ::FREESASA_FAIL.
    
    @param number The residue number will be saved to this string.
    @param line Line from a PDB file.
    @return ::FREESASA_SUCCESS if input is readable, else ::FREESASA_FAIL.
 */
int
freesasa_pdb_get_res_number(char *number,
                            const char* line);

/**
    Get chain label from PDB line.
    
    Extracts the one character chain label (Chain identifier) from an
    `ATOM` or `HETATM` PDB line (i.e. `'A'`, `'B'`, `'C'`, ...) 

    @param line Line from a PDB file.
    @return The chain label. If the line is invalid, the function
    returns '\0'. 
 */
char
freesasa_pdb_get_chain_label(const char* line);

/**
    Get alternate coordinate label from PDB line.
    
    If there is more than one set of coordinates for an atom there
    will be a label 'A', 'B', etc (Alternate location indicator). Else
    the label is ' '. 
    
    @param line Line from a PDB file.
    @return The label. If line is invalid, the function returns
    '\0'.
 */
char
freesasa_pdb_get_alt_coord_label(const char* line);

/**
    Get element symbol from PDB line.

    Writes padded string to argument symbol, i.e. " C", "SE", etc.

    @param symbol The symbol will be written to this string.
    @param line Line from a PDB file.
    @return ::FREESASA_SUCCESS line has more 78 or more characters. ::FREESASA_FAIL if
      line too short.
 */
int
freesasa_pdb_get_symbol(char *symbol,
                        const char* line);

/**
    Is atom Hydrogen?
    
    Checks if an atom from an `ATOM` or `HETATM` PDB line is Hydrogen.

    @param line Line from a PDB file.  
    @return 1 if Hydrogen (or deuterium). 0 otherwise. If line is
    invalid, the function returns ::FREESASA_FAIL.
 */
int
freesasa_pdb_ishydrogen(const char* line);

#endif /* FREESASA_PDB_H */
