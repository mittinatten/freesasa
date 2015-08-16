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

#ifndef FREESASA_CLASSIFY_H
#define FREESASA_CLASSIFY_H

#include <stdio.h>
#include "freesasa.h"

/**
    @file
    @author Simon Mitternacht

    This set of functions maps between different classes of atoms and
    residues. In addition, the function freesasa_classify_radius()
    maps atom type to radius, using the definitions by Ooi et
    al. (1987) for regular protein atoms and by element for other
    atoms.

    The different classes are specified by integers in an interval
    from 0 to n-1. For each class there is a function that returns the
    value of n for that class, e.g. freesasa_classify_nclasses() or
    freesasa_classify_nelements(). Each class also has a function that
    converts the integer to a descriptive string,
    e.g. freesasa_classify_class2str().

    The available classes are (listed as enums)
    - ::freesasa_class (class polar/apolar/...)
    - ::freesasa_residue (residue ALA/CYS/...)
    - ::freesasa_element (element C/N/...)
    - ::freesasa_oons_class (OONS class, plus some specific to FreeSASA):

    The enum for the first class, ::freesasa_class, is in the header
    freesasa.h, because it is used in the main API. The enums for the
    other classes are made pubic in this header, but are only used
    internally for now, and might be expanded as more functionality is
    required.

    The user can also read in a configuration from a file to specify
    their own classes and atomic radii. The functions that do this
    have names that start with freesasa_classify_user_, this
    functionality is still experimental, in that the interface might
    change, but the file format for the config files should be stable
    (except that more features might be added). The file should have
    two sections: `types:` and `atoms:`. The types-section defines
    what types of atoms are available (aliphatic, aromatic, hydroxyl,
    ...), what the radius of that type is and what class a type
    belongs to (polar, apolar, ...). The user is free to define as
    many types and classes as necessary. The atoms-section consists of
    triplets of residue-name, atom-name (as in the corresponding PDB
    entries) and type. A prototype file would be
    
       ~~~
       types:
       C_ALIPHATIC 2.00 apolar
       C_AROMATIC  1.75 apolar
       N 1.55 polar
       
       # this is a comment
       
       atoms:
       ANY N  N             # this is also a comment
       ANY CB C_ALIPHATIC

       ARG CG C_ALIPHATIC

       PRO CB C_AROMATIC
       ~~~

    The residue type `ANY` can be used for atoms that are the same in
    all or most residues (such as backbone atoms). If there is an
    exception for a given amino acid this can be overridden as is
    shown for `PRO CB` in the example.

 */

//! The residue types that are returned by freesasa_classify_residue()
enum freesasa_residue {
    //Regular amino acids
    freesasa_ALA=0, freesasa_ARG, freesasa_ASN, freesasa_ASP,
    freesasa_CYS, freesasa_GLN, freesasa_GLU, freesasa_GLY,
    freesasa_HIS, freesasa_ILE, freesasa_LEU, freesasa_LYS,
    freesasa_MET, freesasa_PHE, freesasa_PRO, freesasa_SER,
    freesasa_THR, freesasa_TRP, freesasa_TYR, freesasa_VAL,
    //some non-standard ones
    freesasa_CSE, freesasa_ASX, freesasa_GLX,
    freesasa_UNK,
    //DNA
    freesasa_DA, freesasa_DC, freesasa_DG, freesasa_DT,
    freesasa_DU, freesasa_DI,
    //RNA
    freesasa_A, freesasa_C, freesasa_G, freesasa_U, freesasa_I, freesasa_T,
    //generic nucleotide
    freesasa_NN
};

//! The element types that are returned by freesasa_classify_element()
enum freesasa_element {
    freesasa_carbon=0, freesasa_oxygen, freesasa_nitrogen,
    freesasa_sulfur, freesasa_phosphorus, freesasa_selenium,
    freesasa_hydrogen, freesasa_element_unknown
};

//! The OONS classes that are returned by freesasa_classify_oons()
enum freesasa_oons_class {
    freesasa_aliphatic_C=0, freesasa_aromatic_C,
    freesasa_carbo_C, freesasa_amide_N,
    freesasa_carbo_O, freesasa_hydroxyl_O,
    freesasa_oons_sulfur,
    freesasa_oons_selenium,
    freesasa_oons_unknown_polar,
    freesasa_oons_unknown
};

/**
    The radius of a given atom type.

    Given an atom this function returns an atomic radius. Regular
    protein atoms are given radii according to OONS (see documentation
    for reference). Atoms in unknown amino acid residues or in nucleic
    acids are given radii based on element (see function
    freesasa_classify_element_radius()). Unknown atom types and
    hydrogens are assigned radius 0.0. 
   
    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return Atom radius in Ångström.
*/
double freesasa_classify_radius(const char *res_name, const char *atom_name);

//////////////////////
// Polar/Apolar/etc //
//////////////////////

/**
    The class of a given atom. 

    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return Atom class (::freesasa_class).
*/
int freesasa_classify_class(const char *res_name, const char *atom_name);

/**
    Class name. 
   
    @param class Class as returned by freesasa_classify_class().
    @return Class name, `NULL` if argument invalid.
 */
const char* freesasa_classify_class2str(int class);

/**
    The number of different classes. 
    
    The total number of available classes, n. The range of
    freesasa_classify_class() is thus [0,n-1]. To facilitate array
    allocation, etc.

    @return Number of classes.
*/
int freesasa_classify_nclasses(void);


/////////////////////////////
// Amino and Nucleic acids //
/////////////////////////////

/**
    Type of residue.

    Returns an integer signifying the type of a residue. There is one
    value for each of the 20 regular amino acids, plus ASX, GLX, XLE,
    CSE and UNK. In addition nucleic acids are treated as separate
    residue types. Return values span from 0 to
    freesasa_classify_nresiduetypes()-1. The standard 20 amino acids
    have indices 0 to 19.
    
    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @return Residue type (::freesasa_residue).
*/
int freesasa_classify_residue(const char *res_name);

/**
    Residue name.

    Returns the name of a residue classified by
    freesasa_classify_residue(). 
    
    @param res Residue type.
    @return Residue name, `NULL` if argument invalid. */
const char* freesasa_classify_residue2str(int res);

/**
    The number of different residue types.

    The total number of available residue types, n. The range of
    freesasa_classify_residue() is thus [0,n-1]. 
    
    @return Number of residue types.
*/
int freesasa_classify_nresiduetypes();


//////////////
// Elements //
//////////////

/**
    Element of a given atom.
    
    Returns an integer signifying the element of an atom. Return
    values span from 0 to freesasa_classify_nelements()-1. 
    
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return The element (::freesasa_element).
*/
int freesasa_classify_element(const char *atom_name);

/**
    Element name.
    
    Returns the name of an element classified by
    freesasa_classify_element(). 
    
    @param element The element.
    @return Element name, `NULL` if argument invalid.
*/
const char* freesasa_classify_element2str(int element);

/**
    The vdW-radius of a given element, according to
    http://www.periodictable.com. 

    @param element The element.  
    @return Radius in Ångström. Unknown
    elements are assigned radius 0.
 */
double freesasa_classify_element_radius(int element);

/**
    The number of element types.

    The total number of available element types, n. The range of
    freesasa_classify_element() is thus [0,n-1].     

    @return Number of element types.
*/
int freesasa_classify_nelements(void);


//////////////////
// OONS classes //
//////////////////

/**
    OONS type of a given atom.

    Returns an integer signifying the type of an atom according to
    OONS (carboxy O, aliphatic C, etc (see preamble to
    classify.h). Return values span from 0 to
    freesasa_classify_noons().  

    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return OONS type (::freesasa_oons_class).
*/
int freesasa_classify_oons(const char *res_name, const char *atom_name);

/**
    OONS type name. 

    @param oons_type The OONS type.
    @return Name of the OONS type, `NULL` if argument invalid.
*/
const char* freesasa_classify_oons2str(int oons_type);

/**
    Class (polar/apolar/etc) for an OONS type. 

    Will not return the class FREESASA_NUCLEIC_ACID, since this is not
    part of the OONS scheme.

    @param oons_type The OONS type.
    @return The class (::freesasa_class).
*/
int freesasa_classify_oons2class(int oons_type);

/**
    Radius of an OONS atom type. 

    @param oons_type The OONS type.
    @return Radius in Ångström. 0.0 Å for hydrogens and unknown
    atoms. */
double freesasa_classify_oons_radius(int oons_type);

/**
    The number of OONS types

    The total number of available OONS types, n. The range of
    freesasa_classify_oobns() is thus [0,n-1].     
  
   @return The number of OONS types.
*/
int freesasa_classify_noons(void);


///////////////////
// Other classes //
///////////////////

/**
    Is residue an amino acid?
    
    Takes the values produced by freesasa_classify_residue() and
    determines if it's an amino acid or not. 

    @param res The residue type.
    @return 1 means amino acid, 0 not, ::FREESASA_FAIL illegal input. 
*/
int freesasa_classify_is_aminoacid(int res);

/**
    Is residue a nucleic acid?
    
    Takes the values produced by freesasa_classify_residue() and
    determines if it's an nucleic acid or not. 

    @param res The residue type.
    @return 1 means nucleic acid, 0 not, ::FREESASA_FAIL illegal input. 
*/
int freesasa_classify_is_nucleicacid(int res);


//////////////
// Validity //
//////////////

/**
    Validate atom description.
    
    Checks if an atom is recognized and can be classified. Prints
    explanatory errors when not recognized.

    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return ::FREESASA_SUCCESS if fully recognized, ::FREESASA_WARN
    if not, and ::FREESASA_FAIL if illegal input. 
*/
int freesasa_classify_validate_atom(const char *res_name, 
                                    const char *atom_name);


#endif
