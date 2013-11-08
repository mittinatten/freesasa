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

#ifndef SASALIB_CLASSIFY_H
#define SASALIB_CLASSIFY_H

#include "sasalib.h"

/** This set of functions maps between different classes of atoms and
    residues. In addition, the function classify_radius(2) maps atom
    type to radius, using the definitions by Ooi et al. (1987) for
    regular protein atoms and by element for other atoms. 

    The different classes are specified by integers in an interval
    from 0 to n-1. For each class there is a function that returns the
    value of n for that class, e.g. classify_nclasses() or
    classify_nelements(). Each class also has a function that converts
    the integer to a descriptive string, e.g. classify_class2str(1).

    The available classes are
       class: polar/apolar/nucleic acid/unknown 
       residue: ALA/ARG/.../ASX/GLX/UNK/A/C/.../DA/DC/.../NN
       element: C/O/N/.../unknown
       oons_type: carbo_C/carbo_O/.../unknown
    
    The enum for the first class, sasalib_class, is visible for the
    user since the list is short and unlikely to change. Also it is
    the one classification where the user is most likely to wish treat
    the different classes individually. The user should not expect the
    numbering, ordering or names of members in the other classes to
    stay constant.
*/

/** Given an atom this function returns an atomic radius. Regular
    protein atoms are given radii according to OONS (see documentation
    for reference). Atoms in unknown amino acid residues or in nucleic
    acids are given radii based on element (see function
    classify_element_radius()). Unknown atom types and hydrogens are
    assigned radius 0.0. */
double classify_radius(const char *res_name, const char *atom_name);

/*********************/
/** Polar/Apolar/etc */
/*********************/

/** Returns the class of a given atom. The return value is one of
    SASALIB_POLAR, SASALIB_APOLAR, SASALIB_NUCLEICACID and
    SASALIB_CLASS_UNKNOWN. */
int classify_class(const char *res_name, const char *atom_name);

/** Takes a class as argument and returns a string that describes that
    class. */
const char* classify_class2str(int class);

/** Returns number of different classes in use. (To facilitate array
    allocation, etc). */
int classify_nclasses();

/****************************/
/** Amino and Nucleic acids */
/****************************/

/** Returns an integer signifying the type of a residue. There is one
    value for each of the 20 regular amino acids, plus ASX, GLX, XLE,
    CSE and UNK. In addition nucleic acids are treated as separate
    residue types. Return values span from 0 to
    classify_nresiduetypes()-1.*/
int classify_residue(const char *res_name);

/** Returns the name of a residue classified by classify_residue().*/
const char* classify_residue2str(int res);

/** Returns the number of possible residue types returned by
    classify_residue(). */
int classify_nresiduetypes();

/**************/
/** Elements **/
/**************/

/** Returns an integer signifying the element of an atom. Return
    values span from 0 to classify_nelements()-1. */
int classify_element(const char *atom_name);

/** Returns a string that describes an element. */
const char* classify_element2str(int element);

/** The vdW-radius of a given element, according to
    www.periodictable.com. */
double classify_element_radius(int element);

/** Returns the number of possible element types returned by
    classify_element(). */
int classify_nelements();

/******************/
/** OONS classes **/
/******************/

/** Returns an integer signifying the type of an atom according to
    OONS (carboxy O, aliphatic C, etc). Return values span from 0 to
    classify_noons(). Hydrogen and Selenium atoms give return value
    sasalib_oons_unknown */
int classify_oons(const char *res_name, const char *atom_name);

/** Returns the name of an OONS-type. */
const char* classify_oons2str(int oons_type);

/** Returns a class (polar/apolar/etc) for an OONS type. Will not
    return the classify_nucleic, since this in not part of the OONS
    scheme. */
int classify_oons2class(int oons_type);

/** Returns the radius of an OONS atom type. Returns 0.0 for hydrogens
    and unknown atoms. */
double classify_oons_radius(int oons_type);

/** Returns the number of available OONS types returned by
    classify_oons(). */
int classify_noons();

/*******************/
/** Other classes **/
/*******************/

/** Takes the values produced by classify_residue() and determines if
    it's an amino acid or not. Return value 1 means amino acid, 0 not,
    -1 illegal input. */
int classify_is_aminoacid(int res);

/** Takes the values produced by classify_residue() and determines if
    it's an nucleic acid or not. Return value 1 is nucleic acid, 0
    not, -1 illegal input. */
int classify_is_nucleicacid(int res);

#endif
