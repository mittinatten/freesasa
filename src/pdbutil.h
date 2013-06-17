#ifndef PDBUTIL_H
#define PDBUTIL_H

#include "vector3.h"

#define PDB_ATOM_NAME_STRL 4
#define PDB_ATOM_RES_NAME_STRL 3
#define PDB_ATOM_RES_NUMBER_STRL 4

/** Extracts the whole atom-name field from an ATOM pdb-line,
    including padding, i.e. a string of PDB_ATOM_NAME_STRL
    characters. For example " CA " for a regular c-alpha */
void pdbutil_get_atom_name(const char *line, char *name);

/** Extracts the whole residue name from an ATOM pdb-line, i.e. a
    string of PDB_ATOM_RES_NAME_STRL characters. For example "ALA" for
    an alanine. */
void pdbutil_get_res_name(const char *line, char *name);

/** Extracts coordinates from an ATOM pdb-line. */
void pdbutil_get_coord(const char *line, vector3 *coord);

/** Extracts residue number as a string from an ATOM pdb-line. String
    format is used because not all residue-numbers are
    numbers. String has length PDB_ATOM_RES_NUMBER_STRL. */
void pdbutil_get_res_number(const char* line, char *number);

/** Extracts the one character chain label from an ATOM pdb-line
    (i.e. A, B, C, ...) */
char pdbutil_get_chain_label(const char* line);

/** Returns 1 if the atom in an ATOM pdb-line is a hydrogen, 0
    otherwise. */
int pdbutil_ishydrogen(const char* line);

#endif
