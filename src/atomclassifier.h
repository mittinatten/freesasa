#ifndef CLASSIFIER_H
#define CLASSIFIER_H

/** this struct can be used to pass information about how to classify atoms
    based on residue and the atom name (in PDB format).
 */

typedef struct {
    int nclasses;
    const char **class2str;
    int (*classify)(const char *res_name, const char *atom_name);
} atomclassifier;

/** This is an example classifier for identifying residues (see
    atomclassifier.c for implementation). In addition to the 20
    standard amino acids it includes ASX, GLX, CSE/SEC and
    PYH/PYL. (CSE and PYH have not been tested, haven't found any
    example PDB files). */
atomclassifier atomclassifier_residue();

#endif
