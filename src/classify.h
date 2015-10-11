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

    The functions in this header are the ones used by
    ::freesasa_default_classifier and
    ::freesasa_residue_classifier. They follow th OONS scheme as far
    as possible and do educated guesses for atoms not included in the 
    original definitions by Ooi et al. (1987). 

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
 */

//! The OONS classes that are returned by freesasa_classify_oons()
enum oons_class {
    oons_aliphatic_C=0, oons_aromatic_C,
    oons_carbo_C, oons_amide_N,
    oons_carbo_O, oons_hydroxyl_O,
    oons_sulfur,
    oons_selenium,
    oons_unknown_polar,
    oons_unknown
};

/**
    The radius of a given atom type.

    Given an atom this function returns an atomic radius. Regular
    protein atoms are given radii according to OONS (see manual
    for reference). Atoms in unknown amino acid residues or in nucleic
    acids are given radii based on element (see function
    freesasa_classify_element_radius()). Unknown atom types and
    hydrogens are assigned radius 0.0. 
   
    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return Atom radius in Ångström.
*/
double
freesasa_classify_radius(const char *res_name,
                         const char *atom_name);

//////////////////////
// Polar/Apolar/etc //
//////////////////////

/**
    The class of a given atom. 

    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @param atom_name The atom name in the format `" CA "`, `" OXT"`, etc.
    @return Atom class (::freesasa_class).
*/
int
freesasa_classify_class(const char *res_name,
                        const char *atom_name);

/**
    Class name. 
   
    @param class Class as returned by freesasa_classify_class().
    @return Class name, `NULL` if argument invalid.
 */
const char*
freesasa_classify_class2str(int class);

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
    `freesasa_classify_nresiduetypes()-1`. The standard 20 amino acids
    have indices 0 to 19.
    
    @param res_name The residue name in the format `"ALA"`, `"PHE"`, etc.
    @return Residue type (::freesasa_residue).
*/
int
freesasa_classify_residue(const char *res_name);

/**
    Residue name.

    Returns the name of a residue classified by
    freesasa_classify_residue(). 
    
    @param res Residue type.
    @return Residue name, `NULL` if argument invalid. */
const char*
freesasa_classify_residue2str(int res);

/**
    The number of different residue types.

    The total number of available residue types, n. The range of
    freesasa_classify_residue() is thus [0,n-1]. 
    
    @return Number of residue types.
*/
int
freesasa_classify_nresiduetypes(void);


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
int
freesasa_classify_element(const char *atom_name);

/**
    Element name.
    
    Returns the name of an element classified by
    freesasa_classify_element(). 
    
    @param element The element.
    @return Element name, `NULL` if argument invalid.
*/
const char*
freesasa_classify_element2str(int element);

/**
    The vdW-radius of a given element, according to
    http://www.periodictable.com. 

    @param element The element.  
    @return Radius in Ångström. Unknown
    elements are assigned radius 0.
 */
double
freesasa_classify_element_radius(int element);

/**
    The number of element types.

    The total number of available element types, n. The range of
    freesasa_classify_element() is thus [0,n-1].     

    @return Number of element types.
*/
int
freesasa_classify_nelements(void);


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
int
freesasa_classify_oons(const char *res_name,
                       const char *atom_name);

/**
    OONS type name. 

    @param oons_type The OONS type.
    @return Name of the OONS type, `NULL` if argument invalid.
*/
const char*
freesasa_classify_oons2str(int oons_type);

/**
    Class (polar/apolar/etc) for an OONS type. 

    Will not return the class FREESASA_NUCLEIC_ACID, since this is not
    part of the OONS scheme.

    @param oons_type The OONS type.
    @return The class (::freesasa_class).
*/
int
freesasa_classify_oons2class(int oons_type);

/**
    Radius of an OONS atom type. 

    @param oons_type The OONS type.
    @return Radius in Ångström. 0.0 Å for hydrogens and unknown
    atoms. */
double
freesasa_classify_oons_radius(int oons_type);

/**
    The number of OONS types

    The total number of available OONS types, n. The range of
    freesasa_classify_oobns() is thus [0,n-1].     
  
   @return The number of OONS types.
*/
int
freesasa_classify_noons(void);


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
int
freesasa_classify_is_aminoacid(int res);

/**
    Is residue a capping group?
    
    Takes the values produced by freesasa_classify_residue() and
    determines if it's a capping group or not. 

    @param res The residue type.
    @return 1 means capping group, 0 not, ::FREESASA_FAIL illegal input. 
 */
int
freesasa_classify_is_capping_group(int res);

/**
    Is residue a nucleic acid?
    
    Takes the values produced by freesasa_classify_residue() and
    determines if it's an nucleic acid or not. 

    @param res The residue type.
    @return 1 means nucleic acid, 0 not, ::FREESASA_FAIL illegal input. 
*/
int
freesasa_classify_is_nucleicacid(int res);


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
int
freesasa_classify_validate_atom(const char *res_name, 
                                const char *atom_name);


#endif /* FREESASA_CLASSIFY_H */
