#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pdbutil.h"
#include "classify.h"

const char *classify_classes[] = {
    "Polar", "Apolar", "Nucleotide", "Unknown"
};

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

enum {ALA=0,ARG,ASN,ASP,
      CYS,GLN,GLU,GLY,
      HIS,ILE,LEU,LYS,
      MET,PHE,PRO,SER,
      THR,TRP,TYR,VAL,
      CSE,ASX,GLX,XLE,
      UNK,
      DA,DC,DG,DT,
      RA,RC,RG,RU,RI,RT,
      NN
};

// the elements seen in PDB Atom entries
const char *element_names[] = {
    "C","O","N","S","P","Se","H","unknown"
};
enum {carbon=0,oxygen,nitrogen,sulphur,phosphorus,selenium,
      hydrogen,element_unknown};

const char *oons_names[] = {
    "aliphatic_C", "aromatic_C",
    "carbo_C", "amide_N", "carbo_O",
    "hydroxyl_O", "sulfur","unknown_polar","unknown"
};
enum {aliphatic_C=0, aromatic_C,
      carbo_C, amide_N, carbo_O, 
      hydroxyl_O, sulfur,
      oons_unknown_polar,
      oons_unknown
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
    return (last-first);
}

double classify_radius(const char *res_name, const char *atom_name)
{
    int res = classify_residue(res_name);
    if (classify_is_aminoacid(res)) {
	return classify_oons_radius(classify_oons(res_name,atom_name));
    } else if (res == UNK) {
	fprintf(stderr,"Warning: residue of type '%s', atom '%s' is unknown."
                " Atom will be classified only according to element.\n",
		res_name,atom_name);
    }
    return classify_element_radius(classify_element(atom_name));
}

int classify_class(const char *res_name, const char *atom_name)
{
    int res = classify_residue(res_name);
    int class;   
    if ((class = classify_oons2class(classify_oons(res_name,atom_name)))
	!= classify_unknown) {
	return class;
    }
    else if (classify_is_nucleicacid(res)) { 
	return classify_nucleicacid;
    }
    return classify_unknown;
}

const char* classify_class2str(int class)
{
    switch (class) {
    case classify_polar: return "Polar";
    case classify_apolar: return "Apolar";
    case classify_nucleicacid: return "Nucleic";
    default: return "Unknown";
    }
}

int classify_nclasses() {
    return classify_unknown+1;
}

int classify_residue(const char *res_name)
{
    char cpy[PDB_ATOM_RES_NAME_STRL+1];
    int len = trim_whitespace(cpy,res_name,strlen(res_name));
    if (len > PDB_ATOM_RES_NAME_STRL) {
	fprintf(stderr,"Warning: Illegal residue name '%s'\n",res_name);
	return UNK;
    }
    for (int i = ALA; i <= NN; ++i) {
        if (! strcmp(res_name,residue_names[i])) return i;
    }

    //for nucleotides we need to deal with whitespace..
    for (int i = DA; i <= NN; ++i) {
	if (! strcmp(cpy,residue_names[i])) return i;
    }

    // warning should perhaps be printed ...
    return UNK;
}

const char* classify_residue2str(int res) {
    if (res < 0 || res > NN) {
	return residue_names[UNK];
    }
    return residue_names[res];
}

int classify_nresiduetypes()
{
    return NN+1;
}

int classify_element(const char *a)
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
    return element_unknown;
}

const char* classify_element2str(int element)
{
    if (element < 0 || element > element_unknown) {
	return element_names[element_unknown];
    }
    return element_names[element];
}

int classify_nelements()
{
    return element_unknown+1;
}

double classify_element_radius(int element)
{
    // based on www.periodictable.com
    assert(element >= 0 && element <= element_unknown);
    switch (element) {
    case carbon: return 1.7;
    case oxygen: return 1.52;
    case nitrogen: return 1.55;
    case sulphur: return 1.8;
    case phosphorus: return 1.8;
    case selenium: return 1.9;
    case hydrogen: return 1.2;
    case element_unknown:
    default:
	return 0.0;
    }
}

/** Functions to deal with the different amino acids in OONS scheme */
int classify_oons_RK(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'N') return amide_N;
    return oons_unknown;
}
int classify_oons_NQDE(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'N') return amide_N;
    if (a[1] == 'O') return carbo_O;
    if (a[1] == 'X') return oons_unknown_polar;
    return oons_unknown;
}
int classify_oons_VIL(const char* a)
{
    if (a[1] == 'C') return aliphatic_C;
    return oons_unknown;
}
int classify_oons_FHPYW(const char* a)
{
    if (a[1] == 'C') return aromatic_C;
    if (a[1] == 'O') return hydroxyl_O;
    if (a[1] == 'N') return amide_N;
    return oons_unknown;
}

int classify_oons_CMST(const char* a) 
{
    if (a[1] == 'C') return aliphatic_C;
    if (a[1] == 'O') return hydroxyl_O;
    if (a[1] == 'S') return sulfur;
    return oons_unknown;
}

int classify_oons_cse(const char* a)
{
    if (a[1] == 'S' && a[2] == 'E') return selenium;
    return oons_unknown;
}

/** Main OONS function */
int classify_oons(const char *res_name, const char *atom_name)
{
    assert(strlen(atom_name) == PDB_ATOM_NAME_STRL);
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);
    
    int res = classify_residue(res_name);

    /* Hydrogens and deuteriums (important to do them here, so they
       can be skipped below */
    if (atom_name[1] == 'H' || atom_name[0] == 'H' ||
	atom_name[1] == 'D' || atom_name[0] == 'D') return hydrogen;

    if (classify_is_aminoacid(res)) {
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
    }
    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    switch (res) {
    case LEU: return classify_oons_VIL(atom_name);
    case SER: return classify_oons_CMST(atom_name);
    case VAL: return classify_oons_VIL(atom_name);
    case GLU: return classify_oons_NQDE(atom_name);
    case LYS: return classify_oons_RK(atom_name);
    case ILE: return classify_oons_VIL(atom_name);
    case THR: return classify_oons_CMST(atom_name);
    case ASP: return classify_oons_NQDE(atom_name);
    case ARG: return classify_oons_RK(atom_name);
    case PRO: return classify_oons_FHPYW(atom_name);
    case ASN: return classify_oons_NQDE(atom_name);
    case PHE: return classify_oons_FHPYW(atom_name);
    case GLN: return classify_oons_NQDE(atom_name);
    case TYR: return classify_oons_FHPYW(atom_name);
    case MET: return classify_oons_CMST(atom_name);
    case HIS: return classify_oons_FHPYW(atom_name);
    case CYS: return classify_oons_CMST(atom_name);
    case TRP: return classify_oons_FHPYW(atom_name);
	// all atoms in Gly and Ala  have already been handled
    case UNK: return oons_unknown;
    case ASX:
    case GLX: return classify_oons_NQDE(atom_name);
    case XLE: return classify_oons_VIL(atom_name);
	// haven't found any PDB files with seleno-cysteine yet,
	// needs testing
	//case CSE: return classify_oons_cse(atom_name);	
    default:
	return oons_unknown;
    }
}

const char* classify_oons2str(int oons_type)
{
    if (oons_type < 0 || oons_type > oons_unknown) {
	return oons_names[oons_unknown];
    }    
    return oons_names[oons_type];
}

int classify_noons()
{
    return oons_unknown+1;
}

int classify_oons2class(int oons_type)
{
    switch (oons_type) {
    case aliphatic_C: return classify_apolar;
    case aromatic_C: return classify_apolar;
    case carbo_C: return classify_apolar;
    case amide_N: return classify_polar;
    case carbo_O: return classify_polar;
    case hydroxyl_O: return classify_polar;
    case sulfur: return classify_polar;
    default: return classify_unknown;
    }
}

double classify_oons_radius(int oons_type)
{
    assert(oons_type >= 0 && oons_type <= oons_unknown);
    switch (oons_type)
    {
    case aliphatic_C: return 2.00;
    case aromatic_C: return 1.75;
    case carbo_C: return 1.55;
    case amide_N: return 1.55;
    case carbo_O: return 1.40;
    case hydroxyl_O: return 1.40;
    case sulfur: return 2.00;
    case oons_unknown:
    default: return 0.0;
    }
}

int classify_is_aminoacid(int res)
{
    if (res >= ALA && res <= XLE) return 1;
    if (res > UNK && res <= NN) return 0;
    return -1;
}

int classify_is_nucleicacid(int res)
{
    if (res >= ALA && res <= UNK) return 0;
    if (res >= DA && res <= NN) return 1;
    return -1;
}
