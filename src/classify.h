#ifndef SASALIB_CLASSIFY_H
#define SASALIB_CLASSIFY_H

/** 4 classes of atoms/chemical groups used*/
typedef enum {
    sasalib_polar=0, sasalib_apolar, 
    sasalib_nucleicacid, sasalib_unknown
} sasalib_class_t;

/** Residue types */
typedef enum {
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
    //RNA
    sasalib_A, sasalib_C, sasalib_G, sasalib_U, sasalib_I, sasalib_T, 
    //generic nucleotide
    sasalib_NN
} sasalib_residue_t;

/** Element types (in ATOM entries, HETATM atom types could be
    anything) */
typedef enum {
    sasalib_carbon=0, sasalib_oxygen, sasalib_nitrogen,
    sasalib_sulfur, sasalib_phosphorus, sasalib_selenium,
    sasalib_hydrogen, sasalib_element_unknown
} sasalib_element_t;

typedef enum {
    sasalib_aliphatic_C=0, sasalib_aromatic_C,
    sasalib_carbo_C, sasalib_amide_N, 
    sasalib_carbo_O, sasalib_hydroxyl_O,
    sasalib_oons_sulfur,
    sasalib_oons_unknown_polar,
    sasalib_oons_unknown
} sasalib_oons_t;

/** Given an atom this function returns an atomic radius. Regular
    protein atoms are given radii according to OONS (see documentation
    for reference). Atoms in unknown amino acid residues or in nucleic
    acids are given radii based on element (see function
    classify_element_radius()). */
double classify_radius(const char *res_name, const char *atom_name);

// Polar/Apolar/etc

/** Returns the class of a given atom. The return value is one of
    sasalib_polar, sasalib_apolar, sasalib_nucleicacid and
    sasalib_unknown. */
sasalib_class_t classify_class(const char *res_name, const char *atom_name);

/** Takes a class as argument and returns a string that describes that
    class. */
const char* classify_class2str(sasalib_class_t class);

/** Returns number of different classes in use. (To facilitate array
    allocation, etc). */
int classify_nclasses();


// Amino and Nucleic acids

/** Returns an integer signifying the type of a residue. There is one
    value for each of the 20 regular amino acids, plus ASX, GLX, XLE,
    CSE and UNK. In addition nucleic acids are treated as separate
    residue types. Return values span from 0 to
    classify_nresiduetypes()-1.*/
sasalib_residue_t classify_residue(const char *res_name);

/** Returns the name of a residue classified by classify_residue().*/
const char* classify_residue2str(sasalib_residue_t res);

/** Returns the number of possible residue types returned by
    classify_residue(). */
int classify_nresiduetypes();


//Elements

/** Returns an integer signifying the element of an atom. Return
    values span from 0 to classify_nelements()-1. */
sasalib_element_t classify_element(const char *atom_name);

/** Returns a string that describes an element. */
const char* classify_element2str(sasalib_element_t element);

/** The vdW-radius of a given element, according to
    www.periodictable.com. */
double classify_element_radius(sasalib_element_t element);

/** Returns the number of possible element types returned by
    classify_element(). */
int classify_nelements();


//OONS classes
 
/** Returns an integer signifying the type of an atom according to
    OONS (carboxy O, aliphatic C, etc). Return values span from 0 to
    classify_noons(). Hydrogen and Selenium atoms give return value
    sasalib_oons_unknown */
sasalib_oons_t classify_oons(const char *res_name, const char *atom_name);

/** Returns the name of an OONS-type. */
const char* classify_oons2str(sasalib_oons_t oons_type);

/** Returns a class (polar/apolar/etc) for an OONS type. Will not
    return the classify_nucleic, since this in not part of the OONS
    scheme. */
sasalib_class_t classify_oons2class(sasalib_oons_t oons_type);

/** Returns the radius of an OONS atom type. Returns 0.0 for hydrogens
    and unknown atoms. */
double classify_oons_radius(sasalib_oons_t oons_type);

/** Returns the number of available OONS types returned by
    classify_oons(). */
int classify_noons();


// Other classes

/** Takes the values produced by classify_residue() and determines if
    it's an amino acid or not. Return value 1 means amino acid, 0 not,
    -1 illegal input. */
int classify_is_aminoacid(sasalib_residue_t res);

/** Takes the values produced by classify_residue() and determines if
    it's an nucleic acid or not. Return value 1 is nucleic acid, 0
    not, -1 illegal input. */
int classify_is_nucleicacid(sasalib_residue_t res);

#endif
