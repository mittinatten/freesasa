#include <string.h>
#include <assert.h>
#include "pdbutil.h"
#include "classify.h"

const char *classify_classes[] = {
    "Polar", "Apolar", "Nucleotide", "Unknown"
};
enum {polar,apolar,nucleotide,class_unknown};

const char *classify_residue_names[] = {
    //amino acids
    "ALA","ARG","ASN","ASP",
    "CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER",
    "TRP","TYR","VAL","THR",
    "CSE","ASX","GLX","UNK"
    //DNA
    "DA","DC","DG","DT",
    //RNA
    "A","C","G","U","I","T"
    //General nuceleotide
    "N"
};
enum {ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,
      PRO,SER,THR,CSE,ASX,GLX,UNK,
      DA,DC,DG,DT,
      A,C,G,U,I,T,
      N
};

// the elements seen in PDB Atom entries
const char *classify_element_names[] = {
    "C","O","N","S","P","Se","H","unknown"
};
enum {carbon,oxygen,nitrogen,sulphur,phosphorus,selenium,
      hydrogen,element_unknown};

const char *classify_oons_names[] = {
    "aliphatic_C", "aromatic_C",
    "carbo_C", "amide_N", "carbo_O",
    "hydroxyl_O", "sulfur","unknown"
};
enum {aliphatic_C, aromatic_C,
      carbo_C, amide_N, carbo_O, 
      hydroxyl_O, sulfur,
      oons_unknown
};

double classify_radius(const char *res_name, const char *atom_name)
{
    if (classify_is_aminoacid(res_name)) {
	return classify_oons_radius(classify_oons(res_name,atom_name));
    } 
    return classify_element_radius(classify_element(atom_name));
}

int classify_class(const char *res_name, const char *atom_name);

const char* classify_class2str(int class);

int classify_residue(const char *res_name, const char *atom_name);

const char* classify_residue2str(int res);

int classify_element(const char *atom_name)
{
    return element_unknown;
}

const char* classify_element2str(int element)
{
    if (element < 0 && element > element_unknown) {
	return classify_element_names[element_unknown];
    }
    return classify_element_names[element];
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

int classify_oons(const char *res_name, const char *atom_name)
{
    return oons_unknown;
}

const char* classify_oons2str(int element);

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
int classify_is_aminoacid(const char *res_name)
{
    return UNK;
}

int classify_is_nucleotide(const char *res_name);

int classify_is_polar(int class);

int classify_is_apolar(int class);
