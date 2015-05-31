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

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "pdb.h"

static inline int pdb_line_check(const char *line,int len) {
    assert(line);
    if (! strncmp(line,"ATOM",4) &&
        ! strncmp(line,"HETATM",6)) {
        return FREESASA_FAIL;
    }
    if (strlen(line) < len) {
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;
}

int freesasa_pdb_get_atom_name(char *name, const char *line)
{
    assert(name);
    assert(line);
    if (pdb_line_check(line,PDB_ATOM_NAME_STRL+12) == FREESASA_FAIL) {
        name[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(name,line+12,PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
    return FREESASA_SUCCESS;
}

int freesasa_pdb_get_res_name(char *name, const char *line)
{
    assert(name);
    assert(line);
    if (pdb_line_check(line,PDB_ATOM_RES_NAME_STRL+17) == FREESASA_FAIL) {
        name[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(name, line+17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
    return FREESASA_SUCCESS;
}
int freesasa_pdb_get_coord(double *xyz, const char *line)
{
    assert(xyz);
    assert(line);
    if (pdb_line_check(line,78) == FREESASA_FAIL) {
        return FREESASA_FAIL;
    }
    sscanf(line+30, "%lf%lf%lf", &xyz[0], &xyz[1], &xyz[2]);
    return FREESASA_SUCCESS;
}
int freesasa_pdb_get_res_number(char *number, const char* line)
{
    assert(number);
    assert(line);
    if (pdb_line_check(line,PDB_ATOM_RES_NUMBER_STRL+22) == FREESASA_FAIL) {
        number[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(number, line+22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
    return FREESASA_SUCCESS;
}
char freesasa_pdb_get_chain_label(const char* line)
{
    assert(line);
    if (pdb_line_check(line,21) == FREESASA_FAIL) return '\0';
    return line[21];
}

char freesasa_pdb_get_alt_coord_label(const char* line)
{
    assert(line);
    if (pdb_line_check(line,16) == FREESASA_FAIL) return '\0';
    return line[16];
}

int freesasa_pdb_ishydrogen(const char* line)
{
    assert(line);
    if (pdb_line_check(line,13) == FREESASA_FAIL) return FREESASA_FAIL;
    //hydrogen
    if (line[12] == 'H' || line[13] == 'H') return 1;
    //hydrogen
    if (line[12] == 'D' || line[13] == 'D') return 1;
    return 0;
}
