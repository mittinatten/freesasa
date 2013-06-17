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

int pdbutil_ishydrogen(const char* line)
{
    if (line[12] == 'H' || line[13] == 'H') return 1;
    return 0;
}
