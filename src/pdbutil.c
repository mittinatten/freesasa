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
#include "pdbutil.h"

void pdbutil_get_atom_name(const char *line, char *name)
{
    strncpy(name,line+12,PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
}

void pdbutil_get_res_name(const char *line, char *name)
{
    strncpy(name, line+17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
}
void pdbutil_get_coord(const char *line, vector3 *coord)
{
    sscanf(line+30, "%lf%lf%lf", &coord->x, &coord->y, &coord->z);
}
void pdbutil_get_res_number(const char* line, char *number)
{
    strncpy(number, line+22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
}
char pdbutil_get_chain_label(const char* line)
{
    return line[21];
}

char pdbutil_get_alt_coord_label(const char* line)
{
    return line[16];
}

int pdbutil_ishydrogen(const char* line)
{
    if (line[12] == 'H' || line[13] == 'H') return 1;
    return 0;
}
