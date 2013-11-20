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

#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "pdb.h"

void sasalib_pdb_get_atom_name(char *name, const char *line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 12+PDB_ATOM_NAME_STRL);
    strncpy(name,line+12,PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
}

void sasalib_pdb_get_res_name(char *name, const char *line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 17+PDB_ATOM_RES_NAME_STRL);
    strncpy(name, line+17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
}
void sasalib_pdb_get_coord(double *xyz, const char *line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 79);
    sscanf(line+30, "%lf%lf%lf", &xyz[0], &xyz[1], &xyz[2]);
}
void sasalib_pdb_get_res_number(char *number, const char* line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 22+PDB_ATOM_RES_NUMBER_STRL);
    strncpy(number, line+22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
}
char sasalib_pdb_get_chain_label(const char* line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 21);
    return line[21];
}

char sasalib_pdb_get_alt_coord_label(const char* line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 16);
    return line[16];
}

int sasalib_pdb_ishydrogen(const char* line)
{
    assert(strlen(line) > 13);
    assert(strncmp(line,"ATOM",4) == 0);
    //hydrogen
    if (line[12] == 'H' || line[13] == 'H') return 1;
    //hydrogen
    if (line[12] == 'D' || line[13] == 'D') return 1;
    return 0;
}
