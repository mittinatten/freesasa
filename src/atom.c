#include <string.h>
#include <assert.h>
#include "atom.h"

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
