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

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "pdb.h"

inline int pdb_line_check(const char *line,int len) {
    if (! strncmp(line,"ATOM",4) && 
	! strncmp(line,"HETATM",6)) {
	return SASALIB_FAIL;
    }
    if (strlen(line) < len) {
	return SASALIB_FAIL;
    }
    return SASALIB_SUCCESS;
}

int sasalib_pdb_get_atom_name(char *name, const char *line)
{
    if (pdb_line_check(line,PDB_ATOM_NAME_STRL+12) == SASALIB_FAIL) {
	name[0] = '\0';
	return SASALIB_FAIL;
    }
    strncpy(name,line+12,PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
    return SASALIB_SUCCESS;
}

int sasalib_pdb_get_res_name(char *name, const char *line)
{
    if (pdb_line_check(line,PDB_ATOM_RES_NAME_STRL+17) == SASALIB_FAIL) {
	name[0] = '\0';
	return SASALIB_FAIL;
    }
    strncpy(name, line+17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
    return SASALIB_SUCCESS;
}
int sasalib_pdb_get_coord(double *xyz, const char *line)
{
    if (pdb_line_check(line,78) == SASALIB_FAIL) {
	return SASALIB_FAIL;
    }
    sscanf(line+30, "%lf%lf%lf", &xyz[0], &xyz[1], &xyz[2]);
    return SASALIB_SUCCESS;
}
int sasalib_pdb_get_res_number(char *number, const char* line)
{
    if (pdb_line_check(line,PDB_ATOM_RES_NUMBER_STRL+22) == SASALIB_FAIL) {
	number[0] = '\0';
	return SASALIB_FAIL;
    }
    strncpy(number, line+22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
    return SASALIB_SUCCESS;
}
char sasalib_pdb_get_chain_label(const char* line)
{
    if (pdb_line_check(line,21) == SASALIB_FAIL) return '\0';
    return line[21];
}

char sasalib_pdb_get_alt_coord_label(const char* line)
{
    if (pdb_line_check(line,16) == SASALIB_FAIL) return '\0';
    return line[16];
}

int sasalib_pdb_ishydrogen(const char* line)
{
    if (pdb_line_check(line,13) == SASALIB_FAIL) return SASALIB_FAIL;
    //hydrogen
    if (line[12] == 'H' || line[13] == 'H') return 1;
    //hydrogen
    if (line[12] == 'D' || line[13] == 'D') return 1;
    return 0;
}
