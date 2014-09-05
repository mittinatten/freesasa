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

#ifndef FREESASA_PDB_H
#define FREESASA_PDB_H

#include "freesasa.h"

#define PDB_ATOM_NAME_STRL 4
#define PDB_ATOM_RES_NAME_STRL 3
#define PDB_ATOM_RES_NUMBER_STRL 4
#define PDB_LINE_STRL 80

/** The following functions all extract info from the PDB lines ATOM
    and HETATM. Valid lines have to begin with either "ATOM" or
    "HETATM" and be sufficiently long to contain the value in
    question. */

/** Extracts the whole atom-name field from an ATOM or HETATM
    pdb-line, including padding, i.e. a string of PDB_ATOM_NAME_STRL
    characters. For example " CA " for a regular c-alpha.  If line is
    invalid, the name will be an empty string and the function returns
    FREESASA_FAIL.  */
int freesasa_pdb_get_atom_name(char *name, const char *line);

/** Extracts the whole residue name from an ATOM or HETATM pdb-line,
    i.e. a string of PDB_ATOM_RES_NAME_STRL characters. For example
    "ALA" for an alanine. If line is invalid, the name will be an
    empty string and the function returns FREESASA_FAIL. */
int freesasa_pdb_get_res_name(char *name, const char *line);

/** Extracts x-, y- and z-coordinates from an ATOM pdb-line. If line
    is invalid, the function returns FREESASA_FAIL and coord will
    remain unchanged. */
int freesasa_pdb_get_coord(double *coord, const char *line);

/** Extracts residue number as a string from an ATOM pdb-line as a
    string. String format is used because not all residue-numbers are
    numbers. The string should have length
    PDB_ATOM_RES_NUMBER_STRL. If line is invalid, number will be an
    empty string and the function returns FREESASA_FAIL.*/
int freesasa_pdb_get_res_number(char *number, const char* line);

/** Extracts the one character chain label from an ATOM pdb-line
    (i.e. A, B, C, ...) If line is invalid, the function returns
    '\0'. */
char freesasa_pdb_get_chain_label(const char* line);

/** If there is more than one set of coordinates for an atom there
    will be a label 'A', 'B', etc. Else the label is ' '. If line is
    invalid, the function returns '\0'. */
char freesasa_pdb_get_alt_coord_label(const char* line);

/** Returns 1 if the atom in an ATOM or HETATM pdb-line is a hydrogen
    (or deuterium), 0 otherwise. If line is invalid, the function returns
    FREESASA_FAIL. */
int freesasa_pdb_ishydrogen(const char* line);

#endif
