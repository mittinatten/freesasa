// why?
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "protein.h"
#include "pdbutil.h"
#include "oons.h"


typedef enum {GLY=0, ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO,
	      SER, THR, CYS, SEC, ASN, GLN, HIS, TYR,
	      ASP, GLU, ARG, LYS, residue_unknown} residue_type;
const int protein_n_residue_types = residue_unknown + 1;

const char *residue_types[] =  
{   "GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO",
    "SER", "THR", "CYS", "SEC", "ASN", "GLN", "HIS", "TYR",
    "ASP", "GLU", "ARG", "LYS", "UNK"};

typedef enum {residue_hydrophobic,residue_polar} residue_class;

struct atom {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char chain_label;
};

struct protein {
    atom *a;
    vector3 *xyz;
    size_t number_atoms;
    size_t number_residues;
    size_t number_chains;
};

const char* residue_type2str(residue_type);

residue_type residue_str2type(const char*);

residue_class residue_type2class(residue_type);

inline residue_class residue_str2class(const char *str)
{
    residue_type2class(residue_str2type(str));
}


protein* protein_init()
{
    protein *p = (protein*) malloc(sizeof(protein));
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = (atom*) malloc(sizeof(atom));
    p->xyz = (vector3*) malloc(sizeof(vector3));
    return p;
}

void protein_free(protein *p)
{
    free(p->a);
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
	    if (pdbutil_ishydrogen(line)) continue;
            vector3 v;
            atom a;
            a.chain_label = pdbutil_get_chain_label(line);
            pdbutil_get_coord(line, &v);
            pdbutil_get_atom_name(line, a.atom_name);
            pdbutil_get_res_name(line, a.res_name);
            pdbutil_get_res_number(line, a.res_number);
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
}

void protein_add_atom(protein *p,
                      const char *atom_name,
                      const char *residue_name,
                      const char *residue_number,
                      char chain_label,
                      double x, double y, double z)
{
    // check input for consistency
    assert(strlen(atom_name) == PDB_ATOM_NAME_STRL);
    assert(strlen(residue_name) == PDB_ATOM_RES_NAME_STRL);
    assert(strlen(residue_number) == PDB_ATOM_RES_NUMBER_STRL);

    // allocate memory, increase number of atoms counter
    protein_alloc_one(p);

    size_t na = p->number_atoms;

    vector3_set(&p->xyz[na-1],x,y,z);

    atom *a = &p->a[na-1];
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    a->chain_label = chain_label;

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

/** get array of radii using custom conversion function, the array r is
    assumed to be of correct size already */
void protein_r(const protein *p, 
	       double *r,
	       double (*atom2radius)(const char *res_name, 
				     const char *atom_name))
{
    for (size_t i = 0; i < p->number_atoms; ++i) {
	r[i] = atom2radius(p->a[i].res_name, p->a[i].atom_name);
    }
}

/** get array of radii using default OONS radii, the array r is
    assumed to of correct size already */
void protein_r_def(const protein *p, double *r)
{
    protein_r(p,r,oons_radius);
}

void protein_update_atom_xyz(protein* p, int number,
			     double x, double y, double z)
{
    assert(number < p->number_atoms);
    vector3_set(&p->xyz[number], x, y, z);
}

void protein_update_atom(protein* p, int number, vector3 *v)
{
    assert(number < p->number_atoms);
    vector3_copy(&p->xyz[number], v);
}

void protein_update_atoms(protein* p, vector3 **v) 
{
    for (int i = 0; i < p->number_atoms; ++i) {
	vector3_copy(&p->xyz[i],v[i]);
    }
}

const vector3* protein_xyz(const protein *p)
{
    return p->xyz;
}

size_t protein_n(const protein *p)
{
    return p->number_atoms;
}

const char* protein_atom_name(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].atom_name;
}

const char* protein_atom_res_name(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_name;
}

const char* protein_atom_res_number(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].res_number;
}

char protein_atom_chain(const protein *p, int i)
{
    assert (i < p->number_atoms);
    return p->a[i].chain_label;
}
