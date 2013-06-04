#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "protein.h"

#define ATOM_NAME_L 4
#define RES_NAME_L 3
#define RES_NUMBER_L 4

void protein_get_atom_name_pdbline(const char *line, char *name)
{
    strncpy(name,line+12,ATOM_NAME_L);
    name[ATOM_NAME_L+1] = '\0';
}

void protein_get_res_name_pdbline(const char *line, char *name)
{
    strncpy(name, line+17, RES_NAME_L);
    name[RES_NAME_L+1] = '\0';
}
void protein_get_coord_pdbline(const char *line, vector3 *coord) 
{
    sscanf(line+30, "%lf%lf%lf", &coord->x, &coord->y, &coord->z);
}
void protein_get_residue_numer(const char* line, char *number)
{
	strncpy(number, line+22, RES_NUMBER_L);
	number[RES_NUMBER_L+1] = '\0';
}
char protein_get_chain_label(const char* line)
{
	return line[21];
}

protein* protein_init()
{
	protein *p = (protein*) malloc(sizeof(protein));
	p->number_atoms = 0;
	p->number_residues = 0;
	p->number_chains = 0;
	return p;
}
void protein_free(protein *p) 
{
	// ...
	free(p->at);
	free(p->ac);
	free(p->coord);
	free(p->chain);
	for (size_t i = 0; i < number_atoms; ++i) {
		free(p->res_name[i]);
		free(p->atom_name[i]);
		free(p->res_number[i]);
	}
	free(p->res_name);
	free(p->atom_name);
	free(p->res_number);
	free(p);
}

protein* protein_init_from_pdb(FILE *pdb_file)
{
	protein* p = protein_init();
	/* two possible implementations, either read file in two passes,
	   first to determine number of atoms, second to read them in,
	   keeping the number of mallocs/reallocs to a minimum.  Or read
	   file once using protein_add_atom() for realloc.  Will begin
	   with second alternative since that has to be implemented
	   anyway. */
	size_t len = 80;
	char line[len];
	while (getline(&line, &len, pdb_file) != -1) {
		if (strncmp("ATOM",line,4)==0) {
			vector3* v;
			char atom_name[5];
			char res_name[5];
			char res_number[5];
			char chain_label = protein_get_chain_pdbline(&line);
			protein_get_coord_pdbline(&line, &v);
			protein_get_atom_name_pdbline(&line, atom_name);
			protein_get_res_name_pdbline(&line, res_name);
			protein_get_res_number_pdbline(&line, res_number);
			protein_add_atom(p,atom_name,res_name,chain_label,
							 res_number, v.x, v.y, v.z);
        }
    }
	//free(line);
	return p;
}

void protein_add_one(protein *p)
{
	int *na = &p->number_atoms;
	++*na;
	p->at = (atom_type*) realloc(p->at,sizeof(atom_type)*(*na));
	p->ac = (atom_class*) realloc(p->ac,sizeof(atom_type)*(*na));
	p->res_name = (char**) realloc(p->res_name,sizeof(char*)*(*na));
	p->atom_name = (char**) realloc(p->atom_name,sizeof(char*)*(*na));
	p->res_number = (char**) realloc(p->atom_name,sizeof(char*)*(*na));
	p->chain = (char*) realloc(p->chain,sizeof(char)*(*na));	
	p->coord = (vector3*) realloc(p->coord,sizeof(vector3)*(*na));
	
	p->res_name[na-1] = (char*) realloc(p->res_name[na-1],
										sizeof(char)*RES_NAME_L);
	p->atom_name[na-1] = (char*) realloc(p->atom_name[na-1],
										 sizeof(char)*ATOM_NAME_L);
	p->res_number[na-1] = (char*) realloc(p->res_number[na-1],
										  sizeof(char)*RES_NUMBER_L);
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
	// possibly check input for consistency
	// we want to keep this fast though..

	// allocate memory, increase number of atoms counter
	protein_add_one(p);

	size_t na = p->number_atoms;
	strcpy(p->atom_name[na-1],atom_name);
	strcpy(p->res_name[na-1],residue_name);
	strcpy(p->res_number[na-1],residue_number);
	p->chain[na-1] = chain_label;
	vector3_set(&p->coord[na-1],x,y,z);

	/* here we assume atoms are ordered sequentially, i.e. chains are
	   not mixed in input, i.e. if two sequential atoms have different
	   residue numbers or chain labels a new residue or new chain is
	   assumed to begin */
	if (p->number_chains == 0) ++p->number_chains; 
	if (na > 2 && chain_label != p->chain[na-2]) ++p->number_chains;

	if (p->number_residues == 0) ++p->number_residues;
	if (na > 2 && strcmp(residue_name,p->chain[na-2])) ++p->number_residues;

	// what else?
	// deal with alternate location indicators...
}
