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

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include "oons.h"
#include "pdbutil.h"

const char *oons_t2s[] = {
    "hydrogen", "aliphatic_C", "aromatic_C",
    "carbo_C", "amide_N", "carbo_O",
    "hydroxyl_O", "sulfur", "selenium",
    "carbon", "oxygen", "nitrogen", "phosphorus",
    "nucleic", "unknown_polar", "unknown"
};

const char *oons_c2s[] = {
    "apolar", "polar", "nucleic", "oons_class_unknown"
};

oons_type oons_name2type (const char *res_name, const char *atom_name);

// regular amino acids
oons_type oons_RK(const char*);
oons_type oons_NQDE(const char*);
oons_type oons_VIL(const char*);
oons_type oons_FHPYW(const char*);
oons_type oons_CMST(const char*);

// nonstandard/uncommon ones
oons_type oons_cse(const char*);

oons_type oons_unk(const char*);

double oons_type2radius(oons_type a)
{
    //OONS Radii [Ooi et al. PNAS 84:3086 (1987)]
    // + elemental vdw-radii for atoms in unknown chemical context
    switch (a)
    {
    case hydrogen: return 0.00; //?? should ideally not happen
    case aliphatic_C: return 2.00;
    case aromatic_C: return 1.75;
    case carbo_C: return 1.55;
    case amide_N: return 1.55;
    case carbo_O: return 1.40;
    case hydroxyl_O: return 1.40;
    case sulfur: return 2.00;
    case carbon:
    case oxygen:
    case nitrogen:
    case phosphorus:
    case type_nucleic: return 1.6; //arbitrary, should be done properly
    case oons_type_unknown:
    default: return 0.0;
    }
}

oons_class oons_type2class(oons_type a)
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
    case type_nucleic: return class_nucleic;
    case hydrogen:
    case oons_type_unknown:
    default: return oons_class_unknown;
    }
}

const char* oons_type2str(oons_type a) 
{
    assert(a <= oons_type_unknown && a >= 0);
    return oons_t2s[a];
}

const char* oons_class2str(oons_class a)
{
    assert(a <= oons_class_unknown && a >= 0);
    return oons_c2s[a];
}

double oons_radius_pdbline(const char *pdb_line)
{
    char res_name[PDB_ATOM_RES_NAME_STRL];
    char atom_name[PDB_ATOM_NAME_STRL];
    pdbutil_get_res_name(pdb_line, res_name);
    pdbutil_get_atom_name(pdb_line, atom_name);
    return oons_radius(res_name, atom_name);
}

double oons_radius(const char *res_name, const char *atom_name) {
    return oons_type2radius(oons_name2type(res_name, atom_name));
}

int oons_classifier(const char *res_name, const char *atom_name)
{
    return oons_type2class(oons_name2type(res_name,atom_name));
}

int oons_typefier(const char *res_name, const char *atom_name)
{
    return oons_name2type(res_name,atom_name);
}

atomclassifier oons_classes()
{
    atomclassifier a;
    a.nclasses = oons_class_unknown+1;
    a.class2str = oons_c2s;
    a.classify = &oons_classifier;
    a.name = "polar_apolar";
    return a;
}

atomclassifier oons_types()
{
    atomclassifier a;
    a.nclasses = oons_type_unknown+1;
    a.class2str = oons_t2s;
    a.classify = &oons_typefier;
    a.name = "OONS_type";
    return a;
}

// support for RNA/DNA should be added at some point..
// pyrrolysine
oons_type oons_name2type (const char *res_name, const char *atom_name)
{
    oons_type type = oons_type_unknown;

    assert(strlen(atom_name) == PDB_ATOM_NAME_STRL);
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);

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

    /* Hydrogens and deuteriums (important to do them here, so they
       can be skipped below */
    if (atom_name[1] == 'H' || atom_name[0] == 'H' ||
	atom_name[1] == 'D' || atom_name[0] == 'D') return hydrogen;

    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    if (! strcmp(res_name, "LEU")) return oons_VIL(atom_name);
    if (! strcmp(res_name, "SER")) return oons_CMST(atom_name);
    if (! strcmp(res_name, "VAL")) return oons_VIL(atom_name);
    if (! strcmp(res_name, "GLU")) return oons_NQDE(atom_name);
    if (! strcmp(res_name, "LYS")) return oons_RK(atom_name);
    if (! strcmp(res_name, "ILE")) return oons_VIL(atom_name);
    if (! strcmp(res_name, "THR")) return oons_CMST(atom_name);
    if (! strcmp(res_name, "ASP")) return oons_NQDE(atom_name);
    if (! strcmp(res_name, "ARG")) return oons_RK(atom_name);
    if (! strcmp(res_name, "PRO")) return oons_FHPYW(atom_name);
    if (! strcmp(res_name, "ASN")) return oons_NQDE(atom_name);
    if (! strcmp(res_name, "PHE")) return oons_FHPYW(atom_name);
    if (! strcmp(res_name, "GLN")) return oons_NQDE(atom_name);
    if (! strcmp(res_name, "TYR")) return oons_FHPYW(atom_name);
    if (! strcmp(res_name, "MET")) return oons_CMST(atom_name);
    if (! strcmp(res_name, "HIS")) return oons_FHPYW(atom_name);
    if (! strcmp(res_name, "CYS")) return oons_CMST(atom_name);
    if (! strcmp(res_name, "TRP")) return oons_FHPYW(atom_name);
    // all atoms in Gly and Ala  have already been handled

    // special cases
    if (! strcmp(res_name, "ASX")) return oons_NQDE(atom_name);
    if (! strcmp(res_name, "GLX")) return oons_NQDE(atom_name);
    if (! strcmp(res_name, "XLE")) return oons_VIL(atom_name);

    // haven't found any PDB files with seleno-cysteine yet, 
    // probably needs work
    //if (! strcmp(res_name, "CSE")) return oons_cse(atom_name);

    //nucleic acids
    if (pdbutil_residuenucleic(res_name)) return class_nucleic;

    if (! strcmp(res_name, "UNK")) {
	fprintf(stderr,"Warning: residue of type 'UNK', atom '%s'. "
                " Atom will be classified only according to element.\n",
		atom_name);
	return oons_unk(atom_name);
    }

    fprintf(stderr,"Warning: combination '%s' '%s' unknown. Atom "
	    "will be assigned radius 0 Å.\n",
	    res_name,atom_name);

    return type;
}

oons_type oons_RK(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'N') return amide_N;
    return oons_type_unknown;
}
oons_type oons_NQDE(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'N') return amide_N;
    if (a[1] == 'O') return carbo_O;
    if (a[1] == 'X') return unknown_polar;
    return oons_type_unknown;
}
oons_type oons_VIL(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    return oons_type_unknown;
}
oons_type oons_FHPYW(const char* a)
{
    if (a[1] == 'C') return aromatic_C;
    if (a[1] == 'O') return hydroxyl_O;
    if (a[1] == 'N') return amide_N;
    return oons_type_unknown;
}

oons_type oons_CMST(const char* a) 
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'O') return hydroxyl_O;
    if (a[1] == 'S') return sulfur;
    return oons_type_unknown;
}

oons_type oons_cse(const char* a)
{
    if (a[1] == 'S' && a[2] == 'E') return selenium;
    return oons_type_unknown;
}

oons_type oons_unk(const char* a)
{
    if ((a[0] == 'C') || (a[0] == ' ' && a[1] == 'C'))
	return carbon;
    if ((a[0] == 'O') || (a[0] == ' ' && a[1] == 'O'))
	return oxygen;
    if ((a[0] == 'N') || (a[0] == ' ' && a[1] == 'N'))
	return nitrogen;
    if ((a[0] == 'S') || (a[0] == ' ' && a[1] == 'S'))
	return sulfur;
    if ((a[0] == 'P') || (a[0] == ' ' && a[1] == 'P'))
	return phosphorus;
    return oons_type_unknown;
}
