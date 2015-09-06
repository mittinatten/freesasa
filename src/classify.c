/*
  Copyright Simon Mitternacht 2013-2015.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include "pdb.h"
#include "classify.h"

extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

// Residue types, make sure this always matches the corresponding enum.
static const char *residue_names[] = {
    //amino acids
    "ALA","ARG","ASN","ASP",
    "CYS","GLN","GLU","GLY",
    "HIS","ILE","LEU","LYS",
    "MET","PHE","PRO","SER",
    "THR","TRP","TYR","VAL",
    "CSE","SEC","ASX","GLX",
    "ACE","NH2",
    "UNK",
    //DNA
    "DA","DC","DG","DT","DU","DI",
    //RNA
    "A","C","G","U","I","T",
    //General nuceleotide
    "N"
        
};


// the elements seen in PDB Atom entries, make sure this always matches the corresponding enum.
static const char *element_names[] = {
    "C","O","N","S","P","Se","H","unknown"
};

// OONS classes, make sure this always matches the corresponding enum.
static const char *oons_names[] = {
    "aliphatic_C", "aromatic_C",
    "carbo_C", "amide_N", "carbo_O",
    "hydroxyl_O", "sulfur","selenium",
    "unknown_polar","unknown"
};

static double default_radius(const char *res_name,
                             const char *atom_name,
                             const freesasa_classifier *c)
{
    return freesasa_classify_radius(res_name,atom_name);
}

static int default_class(const char *res_name,
                         const char *atom_name,
                         const freesasa_classifier *c)
{
    return freesasa_classify_class(res_name,atom_name);
}

const char* default_class2str(int i,
                              const freesasa_classifier *c)
{
    return freesasa_classify_class2str(i);
}

const freesasa_classifier freesasa_default_classifier = {
    .radius = default_radius,
    .sasa_class = default_class,
    .class2str = default_class2str,
    .n_classes = FREESASA_CLASS_UNKNOWN+1,
    .free_config = NULL,
    .config = NULL
};

static int residue(const char *res_name,
                   const char *atom_name,
                   const freesasa_classifier *c)
{
    return freesasa_classify_residue(res_name);
}

static const char* residue2str(int the_residue,
                               const freesasa_classifier *c)
{
    return freesasa_classify_residue2str(the_residue);
}

const freesasa_classifier freesasa_residue_classifier = {
    .radius = default_radius,
    .sasa_class = residue,
    .class2str = residue2str,
    .n_classes = freesasa_NN+1,
    .free_config = NULL,
    .config = NULL
};
  
double freesasa_classify_radius(const char *res_name, const char *atom_name)
{
    if (freesasa_classify_element(atom_name) == freesasa_hydrogen) 
        return freesasa_classify_element_radius(freesasa_hydrogen);
    int res = freesasa_classify_residue(res_name);
    if (freesasa_classify_is_aminoacid(res)) {
        return freesasa_classify_oons_radius(freesasa_classify_oons(res_name,
                                                                    atom_name));
    } else if (res == freesasa_UNK) {
        freesasa_warn("Residue of type '%s', atom '%s' is unknown."
                      " Atom will be classified only according to element.",
                      res_name,atom_name);
    }
    return freesasa_classify_element_radius(freesasa_classify_element(atom_name));
}

int freesasa_classify_class(const char *res_name, const char *atom_name)
{
    int res = freesasa_classify_residue(res_name);
    int class;
    if ((class = freesasa_classify_oons2class(freesasa_classify_oons(res_name,
                                                                     atom_name)))
        != FREESASA_CLASS_UNKNOWN) {
        return class;
    }
    else if (freesasa_classify_is_nucleicacid(res)) {
        return FREESASA_NUCLEICACID;
    }
    return FREESASA_CLASS_UNKNOWN;
}

const char* freesasa_classify_class2str(int class)
{
    switch (class) {
    case FREESASA_POLAR: return "Polar";
    case FREESASA_APOLAR: return "Apolar";
    case FREESASA_NUCLEICACID: return "Nucleic";
    case FREESASA_CLASS_UNKNOWN: return "Unknown";
    default:
        freesasa_warn("Illegal class index '%d' in"
                      "freesasa_classify_class2str(). "
                      "Range is [0,%d].",
                      class,freesasa_classify_nclasses()-1);
        return NULL;
    }
    
}

int freesasa_classify_nclasses() {
    return FREESASA_CLASS_UNKNOWN+1;
}

int freesasa_classify_residue(const char *res_name)
{
    char cpy[PDB_ATOM_RES_NAME_STRL+1];
    if (strlen(res_name) > PDB_ATOM_RES_NAME_STRL) {
        freesasa_warn("%s: unknown residue name '%s' (string too long)",
                      __func__,res_name);
        return freesasa_UNK;
    }
    sscanf(res_name,"%s",cpy);
    for (int i = freesasa_ALA; i <= freesasa_NN; ++i) {
        if (! strcmp(res_name,residue_names[i])) return i;
    }

    //for nucleotides we need to deal with whitespace..
    for (int i = freesasa_DA; i <= freesasa_NN; ++i) {
        if (! strcmp(cpy,residue_names[i])) return i;
    }

    freesasa_warn("%s: residue '%s' unknown.",__func__,res_name);
    return freesasa_UNK;
}

const char* freesasa_classify_residue2str(int res) {
    if (res < 0 || res > freesasa_NN) {
        freesasa_warn("%s: Illegal residue index '%d' passed to "
                      "freesasa_classify_residue2str(1). "
                      "Range is [0,%d]",__func__,res,
                      freesasa_classify_nresiduetypes()-1);
        return NULL;
    }
    return residue_names[res];
}

int freesasa_classify_nresiduetypes()
{
    return freesasa_NN+1;
}

int freesasa_classify_element(const char *atom_name)
{
    //strip whitespace to simplify switch below
    char a[PDB_ATOM_NAME_STRL+1];
    if (strlen(atom_name) > PDB_ATOM_NAME_STRL) {
        freesasa_warn("%s: atom '%s' unknown (string too long).\n",
                      __func__,atom_name);
        return freesasa_element_unknown;
    }
    sscanf(atom_name,"%s",a);
    if (strlen(a) > 0) {
        switch (a[0]) {
        case 'C': return freesasa_carbon;
        case 'O': return freesasa_oxygen;
        case 'N': return freesasa_nitrogen;
        case 'S': return freesasa_sulfur;
        case 'P': return freesasa_phosphorus;
        case 'H': return freesasa_hydrogen;
            // what about Se?
        default:
            freesasa_warn("%s: atom '%s' unknown.\n",__func__,atom_name);
            return freesasa_element_unknown;
        }
    }
    return freesasa_element_unknown;
}

const char* freesasa_classify_element2str(int element)
{
    if (element < 0 || element > freesasa_element_unknown) {
        freesasa_warn("%s: Illegal element index '%d' passed to "
                      "freesasa_classify_element2str(). "
                      "Range is [0,%d]",__func__,element,
                      freesasa_classify_nelements()-1);
        return NULL;
    }
    return element_names[element];
}

int freesasa_classify_nelements()
{
    return freesasa_element_unknown+1;
}

double freesasa_classify_element_radius(int element)
{
    // based on www.periodictable.com
    assert(element >= 0 && element <= freesasa_element_unknown);
    switch (element) {
    case freesasa_carbon: return 1.7;
    case freesasa_oxygen: return 1.52;
    case freesasa_nitrogen: return 1.55;
    case freesasa_sulfur: return 1.8;
    case freesasa_phosphorus: return 1.8;
    case freesasa_selenium: return 1.9;
    case freesasa_hydrogen: return 0.0; // the exception
    case freesasa_element_unknown:
    default:
        return 0.0;
    }
}

/** Functions to deal with the different amino acids in OONS scheme.
    The assumption here is that the backbone plus CB has been handled
    before these functions are called. */
static int classify_oons_RK(const char* a)
{
    if (a[1] == 'C') return freesasa_aliphatic_C;
    if (a[1] == 'N') return freesasa_amide_N;
    return freesasa_oons_unknown;
}

static int classify_oons_ND(const char* a)
{
    if (a[1] == 'C') return freesasa_carbo_C;
    if (a[1] == 'N') return freesasa_amide_N;
    if (a[1] == 'O') return freesasa_carbo_O;
    if (a[1] == 'X') return freesasa_oons_unknown_polar;
    if (a[1] == 'A') return freesasa_oons_unknown_polar;
    return freesasa_oons_unknown;
}
static int classify_oons_QE(const char* a)
{
    if (a[1] == 'N') return freesasa_amide_N;
    if (a[1] == 'O') return freesasa_carbo_O;
    if (a[1] == 'C') {
        if (strcmp(a," CD ") == 0) return freesasa_carbo_C;
        return freesasa_aliphatic_C;
    }
    if (a[1] == 'X') return freesasa_oons_unknown_polar;
    if (a[1] == 'A') return freesasa_oons_unknown_polar;
    return freesasa_oons_unknown;
}

static int classify_oons_VIL(const char* a)
{
    if (a[1] == 'C') return freesasa_aliphatic_C;
    return freesasa_oons_unknown;
}

static int classify_oons_FHPYW(const char* a)
{
    if (a[1] == 'C') return freesasa_aromatic_C;
    if (a[1] == 'O') return freesasa_hydroxyl_O;
    if (a[1] == 'N') return freesasa_amide_N;
    return freesasa_oons_unknown;
}

static int classify_oons_CMST(const char* a)
{
    if (a[1] == 'C') return freesasa_aliphatic_C;
    if (a[1] == 'O') return freesasa_hydroxyl_O;
    if (a[1] == 'S') return freesasa_oons_sulfur;
    return freesasa_oons_unknown;
}

static int classify_oons_cse(const char* a)
{
    if (a[1] == 'S' && a[2] == 'E') return freesasa_oons_selenium;
    return freesasa_oons_unknown;
}

static int classify_oons_nh2(const char* a) 
{
    if (a[1] == 'N' && a[2] == 'H' && a[3] == '2') return freesasa_amide_N;
    return freesasa_oons_unknown;
}

static int classify_oons_ace(const char* a)
{
    if (a[1] == 'C' && a[2] == 'H' && a[3] == '3') return freesasa_aliphatic_C;
    return freesasa_oons_unknown;
}

/** Main OONS function */
int freesasa_classify_oons(const char *res_name, const char *a)
{
    assert(strlen(a) == PDB_ATOM_NAME_STRL);
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);
    
    int res;

    /* Hydrogens and deuteriums (important to do them here, so they
       can be skipped below */
    if (a[1] == 'H' || a[0] == 'H' ||
        a[1] == 'D' || a[0] == 'D') return freesasa_oons_unknown;

    res = freesasa_classify_residue(res_name);

    if (freesasa_classify_is_aminoacid(res)) {
        // backbone
        if (! strcmp(a, " C  ")) return freesasa_carbo_C;
        if (! strcmp(a, " N  ")) return freesasa_amide_N;
        if (! strcmp(a, " CA ")) return freesasa_aliphatic_C;
        if (! strcmp(a, " O  ") ||
            ! strcmp(a, " OXT")) return freesasa_carbo_O;

        // CB is almost always the same
        if (! strcmp(a, " CB ")) {
            if (res == freesasa_PRO) return freesasa_aromatic_C;
            else return freesasa_aliphatic_C;
        }
    }


    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    switch (res) {
    case freesasa_LEU: return classify_oons_VIL(a);
    case freesasa_SER: return classify_oons_CMST(a);
    case freesasa_VAL: return classify_oons_VIL(a);
    case freesasa_GLU: return classify_oons_QE(a);
    case freesasa_LYS: return classify_oons_RK(a);
    case freesasa_ILE: return classify_oons_VIL(a);
    case freesasa_THR: return classify_oons_CMST(a);
    case freesasa_ASP: return classify_oons_ND(a);
    case freesasa_ARG: return classify_oons_RK(a);
    case freesasa_PRO: return classify_oons_FHPYW(a);
    case freesasa_ASN: return classify_oons_ND(a);
    case freesasa_PHE: return classify_oons_FHPYW(a);
    case freesasa_GLN: return classify_oons_QE(a);
    case freesasa_TYR: return classify_oons_FHPYW(a);
    case freesasa_MET: return classify_oons_CMST(a);
    case freesasa_HIS: return classify_oons_FHPYW(a);
    case freesasa_CYS: return classify_oons_CMST(a);
    case freesasa_TRP: return classify_oons_FHPYW(a);
        // all atoms in Gly and Ala  have already been handled
    case freesasa_UNK: return freesasa_oons_unknown;
    case freesasa_ASX: return classify_oons_ND(a);
    case freesasa_GLX: return classify_oons_QE(a);
    case freesasa_CSE: return classify_oons_cse(a);
    case freesasa_SEC: return classify_oons_cse(a);
    case freesasa_ACE: return classify_oons_ace(a);
    case freesasa_NH2: return classify_oons_nh2(a);
    default:
        return freesasa_oons_unknown;
    }
}

const char* freesasa_classify_oons2str(int oons_type)
{
    if (oons_type < 0 || oons_type > freesasa_oons_unknown) {
        freesasa_warn("%s: Illegal OONS type index '%d' passed to "
                      "freesasa_classify_oons2str(1). "
                      "Range is [0,%d]",__func__,oons_type,
                      freesasa_classify_noons()-1);        
        return NULL;
    }
    return oons_names[oons_type];
}

int freesasa_classify_noons()
{
    return freesasa_oons_unknown+1;
}

int freesasa_classify_oons2class(int oons_type)
{
    switch (oons_type) {
    case freesasa_aliphatic_C: return FREESASA_APOLAR;
    case freesasa_aromatic_C: return FREESASA_APOLAR;
    case freesasa_carbo_C: return FREESASA_POLAR;
    case freesasa_amide_N: return FREESASA_POLAR;
    case freesasa_carbo_O: return FREESASA_POLAR;
    case freesasa_hydroxyl_O: return FREESASA_POLAR;
    case freesasa_oons_sulfur: return FREESASA_POLAR;
    case freesasa_oons_selenium: return FREESASA_POLAR;
    case freesasa_oons_unknown_polar: return FREESASA_POLAR;
    default: return FREESASA_CLASS_UNKNOWN;
    }
}

double freesasa_classify_oons_radius(int oons_type)
{
    assert(oons_type >= 0 && oons_type <= freesasa_oons_unknown);
    switch (oons_type)
    {
    case freesasa_aliphatic_C: return 2.00;
    case freesasa_aromatic_C: return 1.75;
    case freesasa_carbo_C: return 1.55;
    case freesasa_amide_N: return 1.55;
    case freesasa_carbo_O: return 1.40;
    case freesasa_hydroxyl_O: return 1.40;
    case freesasa_oons_sulfur: return 2.00;
    case freesasa_oons_selenium: return 1.90;
        //this corresponds to either N or O in ALX and GLX
    case freesasa_oons_unknown_polar: return 1.5;
    case freesasa_oons_unknown:
    default: return 0.0;
    }
}

int freesasa_classify_is_aminoacid(int res)
{
    if (res >= freesasa_ALA && res < freesasa_UNK) return 1;
    if (res >= freesasa_UNK && res <= freesasa_NN) return 0;
    return FREESASA_FAIL;
}

int freesasa_classify_is_nucleicacid(int res)
{
    if (res >= freesasa_ALA && res <= freesasa_UNK) return 0;
    if (res >= freesasa_DA && res <= freesasa_NN) return 1;
    return FREESASA_FAIL;
}

int freesasa_classify_validate_atom(const char *residue_name,
                                    const char *atom_name)
{
    if (strlen(atom_name) != PDB_ATOM_NAME_STRL) {
        return freesasa_fail("%s: Atom name '%s' not valid.",
                             __func__,atom_name);
    }
    if (strlen(residue_name) != PDB_ATOM_RES_NAME_STRL) {
        return freesasa_fail("%s: Residue name '%s' not valid.",
                             __func__,residue_name);
    }
    if (freesasa_classify_class(residue_name,atom_name) != FREESASA_CLASS_UNKNOWN) 
        return FREESASA_SUCCESS;
    if (freesasa_classify_element(atom_name) != freesasa_element_unknown) {
        return FREESASA_SUCCESS;
    } 
    return freesasa_warn("%s: Atom '%s' in '%s' unknown.",
                         __func__,atom_name,residue_name);
}
