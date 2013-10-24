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
#include "pdbutil.h"

void pdbutil_get_atom_name(const char *line, char *name)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 12+PDB_ATOM_NAME_STRL);
    strncpy(name,line+12,PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
}

void pdbutil_get_res_name(const char *line, char *name)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 17+PDB_ATOM_RES_NAME_STRL);
    strncpy(name, line+17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
}
void pdbutil_get_coord(const char *line, vector3 *coord)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 79);
    sscanf(line+30, "%lf%lf%lf", &coord->x, &coord->y, &coord->z);
}
void pdbutil_get_res_number(const char* line, char *number)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 22+PDB_ATOM_RES_NUMBER_STRL);
    strncpy(number, line+22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
}
char pdbutil_get_chain_label(const char* line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 21);
    return line[21];
}

char pdbutil_get_alt_coord_label(const char* line)
{
    assert(strncmp(line,"ATOM",4) == 0);
    assert(strlen(line) > 16);
    return line[16];
}

int pdbutil_ishydrogen(const char* line)
{
    assert(strlen(line) > 13);
    assert(strncmp(line,"ATOM",4) == 0);
    //hydrogen
    if (line[12] == 'H' || line[13] == 'H') return 1;
    //hydrogen
    if (line[12] == 'D' || line[13] == 'D') return 1;
    return 0;
}

int pdbutil_residuenucleic(const char* res_name)
{
    //most common first
    if(! strcmp(res_name, "  A")) return 1;
    if(! strcmp(res_name, "  C")) return 1;
    if(! strcmp(res_name, "  G")) return 1;
    if(! strcmp(res_name, "  T")) return 1;
    if(! strcmp(res_name, "  U")) return 1;
    if(! strcmp(res_name, " DA")) return 1;
    if(! strcmp(res_name, " DC")) return 1;
    if(! strcmp(res_name, " DG")) return 1;
    if(! strcmp(res_name, " DT")) return 1;
    if(! strcmp(res_name, " DU")) return 1;

    //less common
    if(! strcmp(res_name, "  N")) return 1;
    if(! strcmp(res_name, "  I")) return 1;
    if(! strcmp(res_name, " DI")) return 1;

    //alternate formats (not so common?)
    if(! strcmp(res_name, " A ")) return 1;
    if(! strcmp(res_name, " C ")) return 1;
    if(! strcmp(res_name, " G ")) return 1;
    if(! strcmp(res_name, " I ")) return 1;
    if(! strcmp(res_name, " N ")) return 1;
    if(! strcmp(res_name, " T ")) return 1;
    if(! strcmp(res_name, " U ")) return 1;
    return 0;
}

/** Returns 1 if an ATOM pdb-line represents a nucleic acid, 0
    otherwise. */
int pdbutil_isnucleicacid(const char* line) 
{
    char res[PDB_ATOM_RES_NAME_STRL+1];
    pdbutil_get_res_name(line,res);
    return pdbutil_residuenucleic(res);
}
