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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pdb.h"
#include "classify.h"

extern int sasalib_fail(const char *format, ...);
extern int sasalib_warn(const char *format, ...);

/** Residue types */
static const char *residue_names[] = {
    //amino acids
    "ALA","ARG","ASN","ASP",
    "CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    "CSE","ASX","GLX","XLE",
    "UNK",
    //DNA
    "DA","DC","DG","DT","DU","DI",
    //RNA
    "A","C","G","U","I","T",
    //General nuceleotide
    "N"
};
enum {
    //Regular amino acids
    sasalib_ALA=0, sasalib_ARG, sasalib_ASN, sasalib_ASP, 
    sasalib_CYS, sasalib_GLN, sasalib_GLU, sasalib_GLY, 
    sasalib_HIS, sasalib_ILE, sasalib_LEU, sasalib_LYS, 
    sasalib_MET, sasalib_PHE, sasalib_PRO, sasalib_SER, 
    sasalib_THR, sasalib_TRP, sasalib_TYR, sasalib_VAL,
    //some non-standard ones
    sasalib_CSE, sasalib_ASX, sasalib_GLX, sasalib_XLE, 
    sasalib_UNK, 
    //DNA
    sasalib_DA, sasalib_DC, sasalib_DG, sasalib_DT, 
    sasalib_DU, sasalib_DI,  
    //RNA
    sasalib_A, sasalib_C, sasalib_G, sasalib_U, sasalib_I, sasalib_T, 
    //generic nucleotide
    sasalib_NN
};

/** Element types (in ATOM entries, HETATM atom types could be
    anything) */
// the elements seen in PDB Atom entries
static const char *element_names[] = {
    "C","O","N","S","P","Se","H","unknown"
};
enum {
    sasalib_carbon=0, sasalib_oxygen, sasalib_nitrogen,
    sasalib_sulfur, sasalib_phosphorus, sasalib_selenium,
    sasalib_hydrogen, sasalib_element_unknown
};

static const char *oons_names[] = {
    "aliphatic_C", "aromatic_C",
    "carbo_C", "amide_N", "carbo_O",
    "hydroxyl_O", "sulfur","unknown_polar","unknown"
};
enum {
    sasalib_aliphatic_C=0, sasalib_aromatic_C,
    sasalib_carbo_C, sasalib_amide_N, 
    sasalib_carbo_O, sasalib_hydroxyl_O,
    sasalib_oons_sulfur,
    sasalib_oons_unknown_polar,
    sasalib_oons_unknown
};

/** helper function, trims whitespace from beginning and end of string */
size_t sasalib_trim_whitespace(char *target, const char *src, 
			       size_t length)
{
    if (length == 0) { return 0; }
    char *buf = (char*) malloc(length+1);
    char *last = buf+length-1;
    char *first = buf;

    strcpy(buf,src);

    while (*first == ' ') ++first;
    while (*last  == ' ') --last;
    if (first > last) return 0;

    *(last+1) = '\0';

    strcpy(target,first);

    free(buf);
    return (last-first)+1;
}

double sasalib_classify_radius(const char *res_name, const char *atom_name)
{
    int res = sasalib_classify_residue(res_name);
    if (sasalib_classify_is_aminoacid(res)) {
	return sasalib_classify_oons_radius(sasalib_classify_oons(res_name,
								  atom_name));
    } else if (res == sasalib_UNK) {
        sasalib_warn("Residue of type '%s', atom '%s' is unknown."
                     " Atom will be classified only according to element.",
                     res_name,atom_name);
    }
    return sasalib_classify_element_radius(sasalib_classify_element(atom_name));
}

int sasalib_classify_class(const char *res_name, const char *atom_name)
{
    int res = sasalib_classify_residue(res_name);
    int class;   
    if ((class = sasalib_classify_oons2class(sasalib_classify_oons(res_name,
								   atom_name)))
	!= SASALIB_CLASS_UNKNOWN) {
	return class;
    }
    else if (sasalib_classify_is_nucleicacid(res)) { 
	return SASALIB_NUCLEICACID;
    }
    return SASALIB_CLASS_UNKNOWN;
}

const char* sasalib_classify_class2str(int class)
{
    switch (class) {
    case SASALIB_POLAR: return "Polar";
    case SASALIB_APOLAR: return "Apolar";
    case SASALIB_NUCLEICACID: return "Nucleic";
    default: return "Unknown";
    }
}

int sasalib_classify_nclasses() {
    return SASALIB_CLASS_UNKNOWN+1;
}

int sasalib_classify_residue(const char *res_name)
{
    char cpy[PDB_ATOM_RES_NAME_STRL+1];
    int len = sasalib_trim_whitespace(cpy,res_name,strlen(res_name));
    if (len > PDB_ATOM_RES_NAME_STRL) {
        sasalib_warn("unknown residue name '%s'",res_name);
        return sasalib_UNK;
    }
    for (int i = sasalib_ALA; i <= sasalib_NN; ++i) {
        if (! strcmp(res_name,residue_names[i])) return i;
    }

    //for nucleotides we need to deal with whitespace..
    for (int i = sasalib_DA; i <= sasalib_NN; ++i) {
	if (! strcmp(cpy,residue_names[i])) return i;
    }

    sasalib_warn("residue '%s' unknown.",res_name);
    return sasalib_UNK;
}

const char* sasalib_classify_residue2str(int res) {
    if (res < 0 || res > sasalib_NN) {
	return residue_names[sasalib_UNK];
    }
    return residue_names[res];
}

int sasalib_classify_nresiduetypes()
{
    return sasalib_NN+1;
}

int sasalib_classify_element(const char *atom_name)
{
    //strip whitespace to simplify switch below
    char a[PDB_ATOM_NAME_STRL+1];
    int len = sasalib_trim_whitespace(a,atom_name,strlen(atom_name));
    if (len > PDB_ATOM_NAME_STRL) {
	sasalib_warn("atom '%s' unknown.\n",atom_name);
	return sasalib_element_unknown;
    }
    else if (len > 0) {
	switch (a[0]) {
	case 'C': return sasalib_carbon;
	case 'O': return sasalib_oxygen;
	case 'N': return sasalib_nitrogen;
	case 'S': return sasalib_sulfur;
	case 'P': return sasalib_phosphorus;
    // what about Se?
	default: 
	    sasalib_warn("atom '%s' unknown.\n",atom_name);
	    return sasalib_element_unknown;
	}
    }

    return sasalib_element_unknown;
}

const char* sasalib_classify_element2str(int element)
{
    if (element < 0 || element > sasalib_element_unknown) {
	return element_names[sasalib_element_unknown];
    }
    return element_names[element];
}

int sasalib_classify_nelements()
{
    return sasalib_element_unknown+1;
}

double sasalib_classify_element_radius(int element)
{
    // based on www.periodictable.com
    assert(element >= 0 && element <= sasalib_element_unknown);
    switch (element) {
    case sasalib_carbon: return 1.7;
    case sasalib_oxygen: return 1.52;
    case sasalib_nitrogen: return 1.55;
    case sasalib_sulfur: return 1.8;
    case sasalib_phosphorus: return 1.8;
	//case sasalib_selenium: return 1.9;
	//case sasalib_hydrogen: return 1.2;
    case sasalib_element_unknown:
    default:
	return 0.0;
    }
}

/** Functions to deal with the different amino acids in OONS scheme */
static int classify_oons_RK(const char* a)
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    if (a[1] == 'N') return sasalib_amide_N;
    return sasalib_oons_unknown;
}

static int classify_oons_ND(const char* a)
{
    if (a[1] == 'C') return sasalib_carbo_C;
    if (a[1] == 'N') return sasalib_amide_N;
    if (a[1] == 'O') return sasalib_carbo_O;
    if (a[1] == 'X') return sasalib_oons_unknown_polar;
    return sasalib_oons_unknown;
}
static int classify_oons_QE(const char* a)
{
    if (a[1] == 'N') return sasalib_amide_N;
    if (a[1] == 'O') return sasalib_carbo_O;
    if (a[1] == 'X') return sasalib_oons_unknown_polar;
    if (a[1] == 'C') {
	if (strcmp(a," CD ") == 0) return sasalib_carbo_C;
	return sasalib_aliphatic_C;
    }
    return sasalib_oons_unknown;
}

static int classify_oons_VIL(const char* a)
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    return sasalib_oons_unknown;
}

static int classify_oons_FHPYW(const char* a)
{
    if (a[1] == 'C') return sasalib_aromatic_C;
    if (a[1] == 'O') return sasalib_hydroxyl_O;
    if (a[1] == 'N') return sasalib_amide_N;
    return sasalib_oons_unknown;
}

static int classify_oons_CMST(const char* a) 
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    if (a[1] == 'O') return sasalib_hydroxyl_O;
    if (a[1] == 'S') return sasalib_oons_sulfur;
    return sasalib_oons_unknown;
}

static int classify_oons_cse(const char* a)
{
    if (a[1] == 'S' && a[2] == 'E') return sasalib_oons_unknown; //what to do about this
    return sasalib_oons_unknown;
}

/** Main OONS function */
int sasalib_classify_oons(const char *res_name, const char *a)
{
    assert(strlen(a) == PDB_ATOM_NAME_STRL);
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);
    
    /* Hydrogens and deuteriums (important to do them here, so they
       can be skipped below */
    if (a[1] == 'H' || a[0] == 'H' ||
	a[1] == 'D' || a[0] == 'D') return sasalib_oons_unknown;

    int res = sasalib_classify_residue(res_name);

    if (sasalib_classify_is_aminoacid(res)) {
	// backbone
	if (! strcmp(a, " C  ")) return sasalib_carbo_C;
	if (! strcmp(a, " N  ")) return sasalib_amide_N;
	if (! strcmp(a, " CA ")) return sasalib_aliphatic_C;
	if (! strcmp(a, " O  ") ||
	    ! strcmp(a, " OXT")) return sasalib_carbo_O;
	
	// CB is almost always the same
	if (! strcmp(a, " CB ")) {
	    if (res == sasalib_PRO) return sasalib_aromatic_C;
	    else return sasalib_aliphatic_C;
	}
    }


    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    switch (res) {
    case sasalib_LEU: return classify_oons_VIL(a);
    case sasalib_SER: return classify_oons_CMST(a);
    case sasalib_VAL: return classify_oons_VIL(a);
    case sasalib_GLU: return classify_oons_QE(a);
    case sasalib_LYS: return classify_oons_RK(a);
    case sasalib_ILE: return classify_oons_VIL(a);
    case sasalib_THR: return classify_oons_CMST(a);
    case sasalib_ASP: return classify_oons_ND(a);
    case sasalib_ARG: return classify_oons_RK(a);
    case sasalib_PRO: return classify_oons_FHPYW(a);
    case sasalib_ASN: return classify_oons_ND(a);
    case sasalib_PHE: return classify_oons_FHPYW(a);
    case sasalib_GLN: return classify_oons_QE(a);
    case sasalib_TYR: return classify_oons_FHPYW(a);
    case sasalib_MET: return classify_oons_CMST(a);
    case sasalib_HIS: return classify_oons_FHPYW(a);
    case sasalib_CYS: return classify_oons_CMST(a);
    case sasalib_TRP: return classify_oons_FHPYW(a);
	// all atoms in Gly and Ala  have already been handled
    case sasalib_UNK: return sasalib_oons_unknown;
    case sasalib_ASX: return classify_oons_ND(a);
    case sasalib_GLX: return classify_oons_QE(a);
    case sasalib_XLE: return classify_oons_VIL(a);
	// haven't found any PDB files with seleno-cysteine yet,
	// needs testing
    case sasalib_CSE: 
        sasalib_warn("residue type '%s' only has limited support.",	res_name);
	return classify_oons_cse(a);	
    default:
	return sasalib_oons_unknown;
    }
}

const char* sasalib_classify_oons2str(int oons_type)
{
    if (oons_type < 0 || oons_type > sasalib_oons_unknown) {
	return oons_names[sasalib_oons_unknown];
    }    
    return oons_names[oons_type];
}

int sasalib_classify_noons()
{
    return sasalib_oons_unknown+1;
}

int sasalib_classify_oons2class(int oons_type)
{
    switch (oons_type) {
    case sasalib_aliphatic_C: return SASALIB_APOLAR;
    case sasalib_aromatic_C: return SASALIB_APOLAR;
    case sasalib_carbo_C: return SASALIB_APOLAR;
    case sasalib_amide_N: return SASALIB_POLAR;
    case sasalib_carbo_O: return SASALIB_POLAR;
    case sasalib_hydroxyl_O: return SASALIB_POLAR;
    case sasalib_oons_sulfur: return SASALIB_POLAR;
    case sasalib_oons_unknown_polar: return SASALIB_POLAR;
    default: return SASALIB_CLASS_UNKNOWN;
    }
}

double sasalib_classify_oons_radius(int oons_type)
{
    assert(oons_type >= 0 && oons_type <= sasalib_oons_unknown);
    switch (oons_type)
    {
    case sasalib_aliphatic_C: return 2.00;
    case sasalib_aromatic_C: return 1.75;
    case sasalib_carbo_C: return 1.55;
    case sasalib_amide_N: return 1.55;
    case sasalib_carbo_O: return 1.40;
    case sasalib_hydroxyl_O: return 1.40;
    case sasalib_oons_sulfur: return 2.00;
	//this corresponds to either N or O in ALX and GLX
    case sasalib_oons_unknown_polar: return 1.5;
    case sasalib_oons_unknown:
    default: return 0.0;
    }
}

int sasalib_classify_is_aminoacid(int res)
{
    if (res >= sasalib_ALA && res <= sasalib_XLE) return 1;
    if (res >= sasalib_UNK && res <= sasalib_NN) return 0;
    return -1;
}

int sasalib_classify_is_nucleicacid(int res)
{
    if (res >= sasalib_ALA && res <= sasalib_UNK) return 0;
    if (res >= sasalib_DA && res <= sasalib_NN) return 1;
    return -1;
}
