// why?
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "protein.h"

#define ATOM_NAME_STRL 4
#define ATOM_RES_NAME_STRL 3
#define ATOM_RES_NUMBER_STRL 4

typedef enum {
    hydrogen, aliphatic_C, aromatic_C,
    carbo_C, amide_N, carbo_O, 
    hydroxyl_O, sulfur, selenium,
    unknown_polar, atom_type_unknown
} atom_type;

typedef enum {
    apolar, polar, charged, atom_class_unknown
} atom_class;

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
    atom_type at;
    atom_class ac;
    char res_name[ATOM_RES_NAME_STRL+1];
    char res_number[ATOM_RES_NUMBER_STRL+1];
    char atom_name[ATOM_NAME_STRL+1];
    char chain_label;
    vector3 *xyz; //coordinates should be stored elsewhere (for speed)
};

struct protein {
    atom *a;
    vector3 *xyz;
    double *r;
    size_t number_atoms;
    size_t number_residues;
    size_t number_chains;
};

double atom_type2radius(atom_type);

atom_class atom_type2class(atom_type); 

/* Give the padded strings as arguments, i.e. positions 13-16 and
   18-20 of an PDB ATOM entry as atom_name and res_name */
atom_type atom_name2type(const char* res_name, const char* atom_name);

const char* atom_type2str(atom_type);

const char* atom_class2str(atom_class);

const char* residue_type2str(residue_type);

residue_type residue_str2type(const char*);

residue_class residue_type2class(residue_type);

inline residue_class residue_str2class(const char *str)
{
    residue_type2class(residue_str2type(str));
}


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
                      char chain_label,
                      double x, double y, double z)
{
    // check input for consistency
    assert(strlen(atom_name) == ATOM_NAME_STRL);
    assert(strlen(residue_name) == ATOM_RES_NAME_STRL);
    assert(strlen(residue_number) == ATOM_RES_NUMBER_STRL);

    // skip hydrogens and unrecognized atoms
    atom_type at = atom_name2type(residue_name, atom_name);
    if (at == hydrogen || at == atom_type_unknown) return;

    // allocate memory, increase number of atoms counter
    protein_alloc_one(p);

    size_t na = p->number_atoms;

    vector3_set(&p->xyz[na-1],x,y,z);

    atom *a = &p->a[na-1];
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    a->chain_label = chain_label;
    a->at = at;
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

void protein_update_atom_xyz(protein* p, int number,
			     double x, double y, double z)
{
    vector3_set(&p->xyz[number], x, y, z);
}

void protein_update_atom(protein* p, int number, vector3 *v)
{
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

const double* protein_r(const protein *p)
{
    return p->r;
}

size_t protein_n(const protein *p)
{
    return p->number_atoms;
}


/********************
 **** Atoms code ****
 ********************/

// regular amino acids
atom_type atom_RK(const char*);
atom_type atom_NQDE(const char*);
atom_type atom_VIL(const char*);
atom_type atom_FHPYW(const char*);
atom_type atom_CMST(const char*);

// these will probably never be called
atom_type atom_ala(const char*);
atom_type atom_gly(const char*);

// nonstandard/uncommon ones
atom_type atom_sec(const char*);


double atom_type2radius(atom_type a)
{
    //OONS Radii [Ooi et al. PNAS 84:3086 (1987)]
    switch (a)
    {
    case hydrogen: return 0.00; //??
    case aliphatic_C: return 2.00;
    case aromatic_C: return 1.75;
    case carbo_C: return 1.55;
    case amide_N: return 1.55;
    case carbo_O: return 1.40;
    case hydroxyl_O: return 1.40;
    case sulfur: return 2.00;
    case atom_type_unknown:
    default: return 0.0;
    }
}

atom_class atom_type2class(atom_type a)
{
    switch (a)
    {
    case aliphatic_C: return apolar;
    case aromatic_C: return apolar;
    case carbo_C: return apolar;
    case amide_N: return polar;
    case carbo_O: return polar;
    case hydroxyl_O: return polar;
    case sulfur: return polar;
    case hydrogen:
    case atom_type_unknown:
    default: return atom_class_unknown;
    }
}

const char* atom_type2str(atom_type a) 
{
    switch(a)
    {
    case aliphatic_C: return "aliphatic_C";
    case aromatic_C:  return "aromatic_C";
    case carbo_C:     return "carbo_C";
    case amide_N:     return "amide_N";
    case carbo_O:     return "carbo_O";
    case hydroxyl_O:  return "hydroxyl_O";
    case sulfur:      return "sulfur";
    case hydrogen:    return "hydrogen";
    case atom_type_unknown: 
    default: return "unknown";
    }
}

const char* atom_class2str(atom_class a)
{
    switch(a)
    {
    case polar: return "polar";
    case apolar: return "apolar";
    case atom_class_unknown:
    default: return "unknown";	
    }
}


// support for RNA/DNA should be added at some point..
// pyrrolysine
atom_type atom_name2type (const char *res_name, const char *atom_name)
{
    atom_type type = atom_type_unknown;
    assert(strlen(atom_name) == ATOM_NAME_STRL);
    assert(strlen(res_name) == ATOM_RES_NAME_STRL);

    // backbone
    if (! strcmp(atom_name, " C  ")) return carbo_C;
    if (! strcmp(atom_name, " N  ")) return amide_N;
    if (! strcmp(atom_name, " CA ")) return aliphatic_C;
    if (! strcmp(atom_name, " O  ") ||
        ! strcmp(atom_name, " OXT")) return carbo_O;

    // CB is almost always the same
    if (! strcmp(atom_name, " CB ")) {
        if (! strcmp(res_name, "PRO")) return aromatic_C;
        else return aliphatic_C;
    }

    /* Hydrogens (important to do them here, so they can be skipped
       below */
    if (atom_name[1] == 'H' || atom_name[0] == 'H') return hydrogen;

    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    if (! strcmp(res_name, "LEU")) return atom_VIL(atom_name);
    if (! strcmp(res_name, "SER")) return atom_CMST(atom_name);
    if (! strcmp(res_name, "VAL")) return atom_VIL(atom_name);
    if (! strcmp(res_name, "GLU")) return atom_NQDE(atom_name);
    if (! strcmp(res_name, "LYS")) return atom_RK(atom_name);
    if (! strcmp(res_name, "ILE")) return atom_VIL(atom_name);
    if (! strcmp(res_name, "THR")) return atom_CMST(atom_name);
    if (! strcmp(res_name, "ASP")) return atom_NQDE(atom_name);
    if (! strcmp(res_name, "ARG")) return atom_RK(atom_name);
    if (! strcmp(res_name, "PRO")) return atom_FHPYW(atom_name);
    if (! strcmp(res_name, "ASN")) return atom_NQDE(atom_name);
    if (! strcmp(res_name, "PHE")) return atom_FHPYW(atom_name);
    if (! strcmp(res_name, "GLN")) return atom_NQDE(atom_name);
    if (! strcmp(res_name, "TYR")) return atom_FHPYW(atom_name);
    if (! strcmp(res_name, "MET")) return atom_CMST(atom_name);
    if (! strcmp(res_name, "HIS")) return atom_FHPYW(atom_name);
    if (! strcmp(res_name, "CYS")) return atom_CMST(atom_name);
    if (! strcmp(res_name, "TRP")) return atom_FHPYW(atom_name);
    // all atoms in Gly and Ala  have already been handled

    // special cases
    if (! strcmp(res_name, "ASX")) return atom_NQDE(atom_name);
    if (! strcmp(res_name, "GLX")) return atom_NQDE(atom_name);
    if (! strcmp(res_name, "XLE")) return atom_VIL(atom_name);
    // haven't found any PDB files with SEC yet, probably needs work
    if (! strcmp(res_name, "SEC")) return atom_sec(atom_name);

    //need to find PDB file that contains PYL to implement
    //if (! strcmp(res_name, "PYL")) return atom_pyl(atom_name);

    return type;
}

atom_type atom_RK(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'N') return amide_N;
    return atom_type_unknown;
}
atom_type atom_NQDE(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'N') return amide_N;
    if (a[1] == 'O') return carbo_O;
    if (a[1] == 'X') return unknown_polar;
    return atom_type_unknown;
}
atom_type atom_VIL(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    return atom_type_unknown;
}
atom_type atom_FHPYW(const char* a)
{
    if (a[1] == 'C') return aromatic_C;
    if (a[1] == 'O') return hydroxyl_O;
    if (a[1] == 'N') return amide_N;
}

atom_type atom_CMST(const char* a) {
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'O') return hydroxyl_O;
    if (a[1] == 'S') return sulfur;
    return atom_type_unknown;
}

atom_type atom_sec(const char* a)
{
    if (a[1] == 'S' && a[2] == 'E') return selenium;
    return atom_type_unknown;
}

//all atoms handled above
atom_type atom_ala(const char* a)
{
    return atom_type_unknown;
}
atom_type atom_gly(const char* a)
{
    return atom_type_unknown;
}
