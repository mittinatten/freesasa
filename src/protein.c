// why?
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "protein.h"

void protein_get_atom_name_pdbline(const char *line, char *name)
{
    strncpy(name,line+12,ATOM_NAME_STRL);
    name[ATOM_NAME_STRL] = '\0';
}

void protein_get_res_name_pdbline(const char *line, char *name)
{
    strncpy(name, line+17, ATOM_RES_NAME_STRL);
    name[ATOM_RES_NAME_STRL] = '\0';
}
void protein_get_coord_pdbline(const char *line, vector3 *coord)
{
    sscanf(line+30, "%lf%lf%lf", &coord->x, &coord->y, &coord->z);
}
void protein_get_res_number_pdbline(const char* line, char *number)
{
    strncpy(number, line+22, ATOM_RES_NUMBER_STRL);
    number[ATOM_RES_NUMBER_STRL] = '\0';
}
char protein_get_chain_label_pdbline(const char* line)
{
    return line[21];
}

protein* protein_init()
{
    protein *p = (protein*) malloc(sizeof(protein));
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = (atom*) malloc(sizeof(atom));
    p->r = (double*) malloc(sizeof(double));
    p->xyz = (vector3*) malloc(sizeof(vector3));
    return p;
}
void protein_free(protein *p)
{
    free(p->a);
    free(p->r);
    free(p->xyz);
    free(p);
}

protein* protein_init_from_pdb(FILE *pdb_file)
{
    protein *p = protein_init();
    /* two possible implementations, either read file in two passes,
       first to determine number of atoms, second to read them in,
       keeping the number of mallocs/reallocs to a minimum.  Or read
       file once using protein_add_atom() for realloc.  Will begin
       with second alternative since that has to be implemented
       anyway. */
    size_t len = 80;
    char *line = (char*) malloc(sizeof(char)*(len+1));
    while (getline(&line, &len, pdb_file) != -1) {
        if (strncmp("ATOM",line,4)==0) {
            vector3 v;
            atom a;
            a.chain_label = protein_get_chain_label_pdbline(line);
            protein_get_coord_pdbline(line, &v);
            protein_get_atom_name_pdbline(line, a.atom_name);
            protein_get_res_name_pdbline(line, a.res_name);
            protein_get_res_number_pdbline(line, a.res_number);
            protein_add_atom(p,a.atom_name,a.res_name,a.res_number,
                             a.chain_label, v.x, v.y, v.z);
        }
    }
    free(line);
    return p;
}

void protein_alloc_one(protein *p)
{
    size_t na = ++p->number_atoms;
    p->a = (atom*) realloc(p->a,sizeof(atom)*na);
    p->xyz = (vector3*) realloc(p->xyz,sizeof(vector3)*na);
    p->r = (double*) realloc(p->r,sizeof(double)*na);
    p->a[na-1].xyz = &p->xyz[na-1];
}

void protein_add_atom(protein *p,
                      const char *atom_name,
                      const char *residue_name,
                      const char *residue_number,
                      const char chain_label,
                      const double x,
                      const double y,
                      const double z)
{
    // check input for consistency
    assert(strlen(atom_name) == ATOM_NAME_STRL);
    assert(strlen(residue_name) == ATOM_RES_NAME_STRL);
    assert(strlen(residue_number) == ATOM_RES_NUMBER_STRL);

    // allocate memory, increase number of atoms counter
    protein_alloc_one(p);

    size_t na = p->number_atoms;

    vector3_set(&p->xyz[na-1],x,y,z);

    atom *a = &p->a[na-1];
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    a->chain_label = chain_label;
    a->at = atom_name2type(residue_name,atom_name);
    a->ac = atom_type2class(a->at);

    p->r[na-1] = atom_type2radius(a->at);

    /* here we assume atoms are ordered sequentially, i.e. chains are
       not mixed in input: if two sequential atoms have different
       residue numbers or chain labels a new residue or new chain is
       assumed to begin */
    if (p->number_chains == 0)
        ++p->number_chains;
    if (na > 2 && chain_label != p->a[na-2].chain_label)
        ++p->number_chains;

    if (p->number_residues == 0)
        ++p->number_residues;
    if (na > 2 && strcmp(residue_number,p->a[na-2].res_number))
        ++p->number_residues;

    // what else?
    // deal with alternate location indicators...
}
