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
#include "util.h"

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

static double 
default_radius(const char *res_name,
               const char *atom_name,
               const freesasa_classifier *c)
{
    return freesasa_classify_radius(res_name,atom_name);
}

static int
default_class(const char *res_name,
              const char *atom_name,
              const freesasa_classifier *c)
{
    return freesasa_classify_class(res_name,atom_name);
}

const char*
default_class2str(int i,
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

static int
residue(const char *res_name,
        const char *atom_name,
        const freesasa_classifier *c)
{
    return freesasa_classify_residue(res_name);
}

static const char*
residue2str(int the_residue,
            const freesasa_classifier *c)
{
    return freesasa_classify_residue2str(the_residue);
}

const freesasa_classifier freesasa_residue_classifier = {
    .radius = default_radius,
    .sasa_class = residue,
    .class2str = residue2str,
    .n_classes = NN+1,
    .free_config = NULL,
    .config = NULL
};
  
double 
freesasa_classify_radius(const char *res_name,
                         const char *atom_name)
{
    if (freesasa_classify_element(atom_name) == hydrogen) 
        return freesasa_classify_element_radius(hydrogen);
    int res = freesasa_classify_residue(res_name);
    if (freesasa_classify_is_aminoacid(res)) {
        return freesasa_classify_oons_radius(freesasa_classify_oons(res_name,
                                                                    atom_name));
    } else if (res == RES_UNK) {
        freesasa_warn("Residue of type '%s', atom '%s' is unknown."
                      " Atom will be classified only according to element.",
                      res_name,atom_name);
    }
    return freesasa_classify_element_radius(freesasa_classify_element(atom_name));
}

int
freesasa_classify_class(const char *res_name,
                        const char *atom_name)
{
    int res = freesasa_classify_residue(res_name);
    int class;
    if ((class = 
         freesasa_classify_oons2class(freesasa_classify_oons(res_name,
                                                             atom_name)))
        != FREESASA_CLASS_UNKNOWN) {
        return class;
    }
    else if (freesasa_classify_is_nucleicacid(res)) {
        return FREESASA_NUCLEICACID;
    }
    return FREESASA_CLASS_UNKNOWN;
}

const char*
freesasa_classify_class2str(int class)
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

int
freesasa_classify_nclasses(void)
{
    return FREESASA_CLASS_UNKNOWN+1;
}

int
freesasa_classify_residue(const char *res_name)
{
    char cpy[PDB_ATOM_RES_NAME_STRL+1];
    if (strlen(res_name) > PDB_ATOM_RES_NAME_STRL) {
        freesasa_warn("%s: unknown residue name '%s' (string too long)",
                      __func__,res_name);
        return RES_UNK;
    }
    sscanf(res_name,"%s",cpy);
    for (int i = ALA; i <= NN; ++i) {
        if (! strcmp(res_name,residue_names[i])) return i;
    }

    //for nucleotides we need to deal with whitespace..
    for (int i = DA; i <= NN; ++i) {
        if (! strcmp(cpy,residue_names[i])) return i;
    }

    freesasa_warn("%s: residue '%s' unknown.",__func__,res_name);
    return RES_UNK;
}

const char*
freesasa_classify_residue2str(int res) 
{
    if (res < 0 || res > NN) {
        freesasa_warn("%s: Illegal residue index '%d' passed to "
                      "freesasa_classify_residue2str(1). "
                      "Range is [0,%d]",__func__,res,
                      freesasa_classify_nresiduetypes()-1);
        return NULL;
    }
    return residue_names[res];
}

int
freesasa_classify_nresiduetypes(void)
{
    return NN+1;
}

int
freesasa_classify_element(const char *atom_name)
{
    //strip whitespace to simplify switch below
    char a[PDB_ATOM_NAME_STRL+1];
    if (strlen(atom_name) > PDB_ATOM_NAME_STRL) {
        freesasa_warn("%s: atom '%s' unknown (string too long).\n",
                      __func__,atom_name);
        return element_unknown;
    }
    sscanf(atom_name,"%s",a);
    if (strlen(a) > 0) {
        switch (a[0]) {
        case 'C': return carbon;
        case 'O': return oxygen;
        case 'N': return nitrogen;
        case 'S': return sulfur;
        case 'P': return phosphorus;
        case 'H': return hydrogen;
            // what about Se?
        default:
            freesasa_warn("%s: atom '%s' unknown.\n",__func__,atom_name);
            return element_unknown;
        }
    }
    return element_unknown;
}

const char*
freesasa_classify_element2str(int element)
{
    if (element < 0 || element > element_unknown) {
        freesasa_warn("%s: Illegal element index '%d' passed to "
                      "freesasa_classify_element2str(). "
                      "Range is [0,%d]",__func__,element,
                      freesasa_classify_nelements()-1);
        return NULL;
    }
    return element_names[element];
}

int
freesasa_classify_nelements()
{
    return element_unknown+1;
}

double
freesasa_classify_element_radius(int element)
{
    // based on www.periodictable.com
    assert(element >= 0 && element <= element_unknown);
    switch (element) {
    case carbon: return 1.7;
    case oxygen: return 1.52;
    case nitrogen: return 1.55;
    case sulfur: return 1.8;
    case phosphorus: return 1.8;
    case selenium: return 1.9;
    case hydrogen: return 0.0; // the exception
    case element_unknown:
    default:
        return 0.0;
    }
}

/** Functions to deal with the different amino acids in OONS scheme.
    The assumption here is that the backbone plus CB has been handled
    before these functions are called. */
static int
oons_RK(const char* a)
{
    if (a[1] == 'C') return oons_aliphatic_C;
    if (a[1] == 'N') return oons_amide_N;
    return oons_unknown;
}

static int
oons_ND(const char* a)
{
    if (a[1] == 'C') return oons_carbo_C;
    if (a[1] == 'N') return oons_amide_N;
    if (a[1] == 'O') return oons_carbo_O;
    if (a[1] == 'X') return oons_unknown_polar;
    if (a[1] == 'A') return oons_unknown_polar;
    return oons_unknown;
}

static int
oons_QE(const char* a)
{
    if (a[1] == 'N') return oons_amide_N;
    if (a[1] == 'O') return oons_carbo_O;
    if (a[1] == 'C') {
        if (strcmp(a," CD ") == 0) return oons_carbo_C;
        return oons_aliphatic_C;
    }
    if (a[1] == 'X') return oons_unknown_polar;
    if (a[1] == 'A') return oons_unknown_polar;
    return oons_unknown;
}

static int
oons_VIL(const char* a)
{
    if (a[1] == 'C') return oons_aliphatic_C;
    return oons_unknown;
}

static int
oons_FHPYW(const char* a)
{
    if (a[1] == 'C') return oons_aromatic_C;
    if (a[1] == 'O') return oons_hydroxyl_O;
    if (a[1] == 'N') return oons_amide_N;
    return oons_unknown;
}

static int
oons_CMST(const char* a)
{
    if (a[1] == 'C') return oons_aliphatic_C;
    if (a[1] == 'O') return oons_hydroxyl_O;
    if (a[1] == 'S') return oons_sulfur;
    return oons_unknown;
}

static int
oons_cse(const char* a)
{
    if (a[0] == 'S' && a[1] == 'E') return oons_selenium;
    return oons_unknown;
}

static int
oons_nh2(const char* a) 
{
    if (a[1] == 'N' && a[2] == 'H' && a[3] == '2') return oons_amide_N;
    return oons_unknown;
}

static int
oons_ace(const char* a)
{
    if (a[1] == 'C' && a[2] == 'H' && a[3] == '3') return oons_aliphatic_C;
    return oons_unknown;
}

/** Main OONS function */
int
freesasa_classify_oons(const char *res_name,
                       const char *a)
{
    assert(strlen(a) == PDB_ATOM_NAME_STRL);
    assert(strlen(res_name) == PDB_ATOM_RES_NAME_STRL);
    
    int res;

    /* Hydrogens and deuteriums (important to do them here, so they
       can be skipped below */
    if (a[1] == 'H' || a[0] == 'H' ||
        a[1] == 'D' || a[0] == 'D') return oons_unknown;

    res = freesasa_classify_residue(res_name);

    if (freesasa_classify_is_aminoacid(res)) {
        // backbone
        if (! strcmp(a, " C  ")) return oons_carbo_C;
        if (! strcmp(a, " N  ")) return oons_amide_N;
        if (! strcmp(a, " CA ")) return oons_aliphatic_C;
        if (! strcmp(a, " O  ") ||
            ! strcmp(a, " OXT")) return oons_carbo_O;

        // CB is almost always the same
        if (! strcmp(a, " CB ")) {
            if (res == PRO) return oons_aromatic_C;
            else return oons_aliphatic_C;
        }
    }


    /* Amino acids are sorted by frequency of occurence for
       optimization (probably has minimal effect, but easy to do) */
    switch (res) {
    case LEU: return oons_VIL(a);
    case SER: return oons_CMST(a);
    case VAL: return oons_VIL(a);
    case GLU: return oons_QE(a);
    case LYS: return oons_RK(a);
    case ILE: return oons_VIL(a);
    case THR: return oons_CMST(a);
    case ASP: return oons_ND(a);
    case ARG: return oons_RK(a);
    case PRO: return oons_FHPYW(a);
    case ASN: return oons_ND(a);
    case PHE: return oons_FHPYW(a);
    case GLN: return oons_QE(a);
    case TYR: return oons_FHPYW(a);
    case MET: return oons_CMST(a);
    case HIS: return oons_FHPYW(a);
    case CYS: return oons_CMST(a);
    case TRP: return oons_FHPYW(a);
        // all atoms in Gly and Ala  have already been handled
    case RES_UNK: return oons_unknown;
    case ASX: return oons_ND(a);
    case GLX: return oons_QE(a);
    case CSE: return oons_cse(a);
    case SEC: return oons_cse(a);
    case ACE: return oons_ace(a);
    case NH2: return oons_nh2(a);
    default:
        return oons_unknown;
    }
}

const char*
freesasa_classify_oons2str(int oons_type)
{
    if (oons_type < 0 || oons_type > oons_unknown) {
        freesasa_warn("%s: Illegal OONS type index '%d' passed to "
                      "freesasa_classify_oons2str(1). "
                      "Range is [0,%d]",__func__,oons_type,
                      freesasa_classify_noons()-1);        
        return NULL;
    }
    return oons_names[oons_type];
}

int
freesasa_classify_noons()
{
    return oons_unknown+1;
}

int
freesasa_classify_oons2class(int oons_type)
{
    switch (oons_type) {
    case oons_aliphatic_C: return FREESASA_APOLAR;
    case oons_aromatic_C: return FREESASA_APOLAR;
    case oons_carbo_C: return FREESASA_POLAR;
    case oons_amide_N: return FREESASA_POLAR;
    case oons_carbo_O: return FREESASA_POLAR;
    case oons_hydroxyl_O: return FREESASA_POLAR;
    case oons_sulfur: return FREESASA_POLAR;
    case oons_selenium: return FREESASA_POLAR;
    case oons_unknown_polar: return FREESASA_POLAR;
    default: return FREESASA_CLASS_UNKNOWN;
    }
}

double
freesasa_classify_oons_radius(int oons_type)
{
    assert(oons_type >= 0 && oons_type <= oons_unknown);
    switch (oons_type)
    {
    case oons_aliphatic_C: return 2.00;
    case oons_aromatic_C: return 1.75;
    case oons_carbo_C: return 1.55;
    case oons_amide_N: return 1.55;
    case oons_carbo_O: return 1.40;
    case oons_hydroxyl_O: return 1.40;
    case oons_sulfur: return 2.00;
    case oons_selenium: return 1.90;
        //this corresponds to either N or O in ALX and GLX
    case oons_unknown_polar: return 1.5;
    case oons_unknown:
    default: return 0.0;
    }
}

int
freesasa_classify_is_aminoacid(int res)
{
    if (res >= ALA && res < RES_UNK) return 1;
    if (res >= RES_UNK && res <= NN) return 0;
    return FREESASA_FAIL;
}

int
freesasa_classify_is_nucleicacid(int res)
{
    if (res >= ALA && res <= RES_UNK) return 0;
    if (res >= DA && res <= NN) return 1;
    return FREESASA_FAIL;
}

int
freesasa_classify_validate_atom(const char *residue_name,
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
    if (freesasa_classify_element(atom_name) != element_unknown) {
        return FREESASA_SUCCESS;
    } 
    return freesasa_warn("%s: Atom '%s' in '%s' unknown.",
                         __func__,atom_name,residue_name);
}
