#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pdbutil.h"
#include "classify.h"



const char *residue_names[] = {
    //amino acids
    "ALA","ARG","ASN","ASP",
    "CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    "CSE","ASX","GLX","XLE",
    "UNK",
    //DNA
    "DA","DC","DG","DT",
    //RNA
    "A","C","G","U","I","T",
    //General nuceleotide
    "N"
};

// the elements seen in PDB Atom entries
const char *element_names[] = {
    "C","O","N","S","P","Se","H","unknown"
};

const char *oons_names[] = {
    "aliphatic_C", "aromatic_C",
    "carbo_C", "amide_N", "carbo_O",
    "hydroxyl_O", "sulfur","unknown_polar","unknown"
};

/** helper function, trims whitespace from beginning and end of string */
size_t trim_whitespace(char *target, const char *src, size_t length)
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

double classify_radius(const char *res_name, const char *atom_name)
{
    int res = classify_residue(res_name);
    if (classify_is_aminoacid(res)) {
	return classify_oons_radius(classify_oons(res_name,atom_name));
    } else if (res == sasalib_UNK) {
	fprintf(stderr,"Warning: residue of type '%s', atom '%s' is unknown."
                " Atom will be classified only according to element.\n",
		res_name,atom_name);
    }
    return classify_element_radius(classify_element(atom_name));
}

sasalib_class_t classify_class(const char *res_name, const char *atom_name)
{
    int res = classify_residue(res_name);
    int class;   
    if ((class = classify_oons2class(classify_oons(res_name,atom_name)))
	!= sasalib_unknown) {
	return class;
    }
    else if (classify_is_nucleicacid(res)) { 
	return sasalib_nucleicacid;
    }
    return sasalib_unknown;
}

const char* classify_class2str(sasalib_class_t class)
{
    switch (class) {
    case sasalib_polar: return "Polar";
    case sasalib_apolar: return "Apolar";
    case sasalib_nucleicacid: return "Nucleic";
    default: return "Unknown";
    }
}

int classify_nclasses() {
    return sasalib_unknown+1;
}

sasalib_residue_t classify_residue(const char *res_name)
{
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);
    char cpy[PDB_ATOM_RES_NAME_STRL+1];
    int len = trim_whitespace(cpy,res_name,strlen(res_name));
    if (len > PDB_ATOM_RES_NAME_STRL) {
	fprintf(stderr,"Warning: Illegal residue name '%s'\n",res_name);
	return sasalib_UNK;
    }
    for (int i = sasalib_ALA; i <= sasalib_NN; ++i) {
        if (! strcmp(res_name,residue_names[i])) return i;
    }

    //for nucleotides we need to deal with whitespace..
    for (int i = sasalib_DA; i <= sasalib_NN; ++i) {
	if (! strcmp(cpy,residue_names[i])) return i;
    }

    fprintf(stderr,"Warning: Residue '%s' unknown.\n",res_name);
    return sasalib_UNK;
}

const char* classify_residue2str(sasalib_residue_t res) {
    if (res < 0 || res > sasalib_NN) {
	return residue_names[sasalib_UNK];
    }
    return residue_names[res];
}

int classify_nresiduetypes()
{
    return sasalib_NN+1;
}

sasalib_element_t classify_element(const char *atom_name)
{
    assert(strlen(atom_name) == PDB_ATOM_NAME_STRL);

    //strip whitespace to simplify switch below
    char a[PDB_ATOM_NAME_STRL+1];
    int len = trim_whitespace(a,atom_name,strlen(atom_name));

    if (len > 0) {
	switch (a[0]) {
	case 'C': return sasalib_carbon;
	case 'O': return sasalib_oxygen;
	case 'N': return sasalib_nitrogen;
	case 'S': return sasalib_sulfur;
	case 'P': return sasalib_phosphorus;
    // what about Se?
	default: 
	    fprintf(stderr,"Warning: Atom '%s' unknown.\n",atom_name);
	    return sasalib_element_unknown;
	}
    }

    return sasalib_element_unknown;
}

const char* classify_element2str(sasalib_element_t element)
{
    if (element < 0 || element > sasalib_element_unknown) {
	return element_names[sasalib_element_unknown];
    }
    return element_names[element];
}

int classify_nelements()
{
    return sasalib_element_unknown+1;
}

double classify_element_radius(sasalib_element_t element)
{
    // based on www.periodictable.com
    assert(element >= 0 && element <= sasalib_element_unknown);
    switch (element) {
    case sasalib_carbon: return 1.7;
    case sasalib_oxygen: return 1.52;
    case sasalib_nitrogen: return 1.55;
    case sasalib_sulfur: return 1.8;
    case sasalib_phosphorus: return 1.8;
    case sasalib_selenium: return 1.9;
    case sasalib_hydrogen: return 1.2;
    case sasalib_element_unknown:
    default:
	return 0.0;
    }
}

/** Functions to deal with the different amino acids in OONS scheme */
int classify_oons_RK(const char* a)
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    if (a[1] == 'N') return sasalib_amide_N;
    return sasalib_oons_unknown;
}
int classify_oons_NQDE(const char* a)
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    if (a[1] == 'N') return sasalib_amide_N;
    if (a[1] == 'O') return sasalib_carbo_O;
    if (a[1] == 'X') return sasalib_oons_unknown_polar;
    return sasalib_oons_unknown;
}
int classify_oons_VIL(const char* a)
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    return sasalib_oons_unknown;
}
int classify_oons_FHPYW(const char* a)
{
    if (a[1] == 'C') return sasalib_aromatic_C;
    if (a[1] == 'O') return sasalib_hydroxyl_O;
    if (a[1] == 'N') return sasalib_amide_N;
    return sasalib_oons_unknown;
}

int classify_oons_CMST(const char* a) 
{
    if (a[1] == 'C') return sasalib_aliphatic_C;
    if (a[1] == 'O') return sasalib_hydroxyl_O;
    if (a[1] == 'S') return sasalib_oons_sulfur;
    return sasalib_oons_unknown;
}

int classify_oons_cse(const char* a)
{
    if (a[1] == 'S' && a[2] == 'E') return sasalib_oons_unknown; //what to do about this
    return sasalib_oons_unknown;
}

/** Main OONS function */
sasalib_oons_t classify_oons(const char *res_name, const char *a)
{
    assert(strlen(a) == PDB_ATOM_NAME_STRL);
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);
    
    sasalib_residue_t res = classify_residue(res_name);

    /* Hydrogens and deuteriums (important to do them here, so they
       can be skipped below */
    if (a[1] == 'H' || a[0] == 'H' ||
	a[1] == 'D' || a[0] == 'D') return sasalib_oons_unknown;

    if (classify_is_aminoacid(res)) {
	// backbone
	if (! strcmp(a, " C  ")) return sasalib_carbo_C;
	if (! strcmp(a, " N  ")) return sasalib_amide_N;
	if (! strcmp(a, " CA ")) return sasalib_aliphatic_C;
	if (! strcmp(a, " O  ") ||
	    ! strcmp(a, " OXT")) return sasalib_carbo_O;
	
	// CB is almost always the same
	if (! strcmp(a, " CB ")) {
	    if (! strcmp(res_name, "PRO")) return sasalib_aromatic_C;
	    else return sasalib_aliphatic_C;
	}
    }
    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    switch (res) {
    case sasalib_LEU: return classify_oons_VIL(a);
    case sasalib_SER: return classify_oons_CMST(a);
    case sasalib_VAL: return classify_oons_VIL(a);
    case sasalib_GLU: return classify_oons_NQDE(a);
    case sasalib_LYS: return classify_oons_RK(a);
    case sasalib_ILE: return classify_oons_VIL(a);
    case sasalib_THR: return classify_oons_CMST(a);
    case sasalib_ASP: return classify_oons_NQDE(a);
    case sasalib_ARG: return classify_oons_RK(a);
    case sasalib_PRO: return classify_oons_FHPYW(a);
    case sasalib_ASN: return classify_oons_NQDE(a);
    case sasalib_PHE: return classify_oons_FHPYW(a);
    case sasalib_GLN: return classify_oons_NQDE(a);
    case sasalib_TYR: return classify_oons_FHPYW(a);
    case sasalib_MET: return classify_oons_CMST(a);
    case sasalib_HIS: return classify_oons_FHPYW(a);
    case sasalib_CYS: return classify_oons_CMST(a);
    case sasalib_TRP: return classify_oons_FHPYW(a);
	// all atoms in Gly and Ala  have already been handled
    case sasalib_UNK: return sasalib_oons_unknown;
    case sasalib_ASX:
    case sasalib_GLX: return classify_oons_NQDE(a);
    case sasalib_XLE: return classify_oons_VIL(a);
	// haven't found any PDB files with seleno-cysteine yet,
	// needs testing
	//case sasalib_CSE: return classify_oons_cse(a);	
    default:
	return sasalib_oons_unknown;
    }
}

const char* classify_oons2str(sasalib_oons_t oons_type)
{
    if (oons_type < 0 || oons_type > sasalib_oons_unknown) {
	return oons_names[sasalib_oons_unknown];
    }    
    return oons_names[oons_type];
}

int classify_noons()
{
    return sasalib_oons_unknown+1;
}

sasalib_class_t classify_oons2class(sasalib_oons_t oons_type)
{
    switch (oons_type) {
    case sasalib_aliphatic_C: return sasalib_apolar;
    case sasalib_aromatic_C: return sasalib_apolar;
    case sasalib_carbo_C: return sasalib_apolar;
    case sasalib_amide_N: return sasalib_polar;
    case sasalib_carbo_O: return sasalib_polar;
    case sasalib_hydroxyl_O: return sasalib_polar;
    case sasalib_oons_sulfur: return sasalib_polar;
    case sasalib_oons_unknown_polar: return sasalib_polar;
    default: return sasalib_unknown;
    }
}

double classify_oons_radius(sasalib_oons_t oons_type)
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
    case sasalib_oons_unknown:
    default: return 0.0;
    }
}

int classify_is_aminoacid(sasalib_residue_t res)
{
    if (res >= sasalib_ALA && res <= sasalib_XLE) return 1;
    if (res >= sasalib_UNK && res <= sasalib_NN) return 0;
    return -1;
}

int classify_is_nucleicacid(sasalib_residue_t res)
{
    if (res >= sasalib_ALA && res <= sasalib_UNK) return 0;
    if (res >= sasalib_DA && res <= sasalib_NN) return 1;
    return -1;
}
