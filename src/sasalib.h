/*
  Copyright Simon Mitternacht 2013-2014.

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

#ifndef SASALIB_H
#define SASALIB_H

/** sasalib.h provides functions to init and perform a SASA
    calculation using Sasalib with standard atom classifications. The
    user may optionally select algorithm and provide parameters. 
    
    * Coordinates *

    If users wish to supply their own coordinates and radii, these are
    accepted as arrays of doubles. The coordinate-array should have
    size 3*n with coordinates in the order x1,y1,z1,x2,y2,z2,... 

    * Error-reporting * 

    All errors are written to stderr and are prefixed with the string
    'sasalib'. There are two error return values SASALIB_WARN and
    SASALIB_FAIL (see documentation of each function to see when these
    are used). The former refers to minor errors where a calculations
    can still proceed, but details might not be as intended
    (e.g. default parameters are used instead of those provided by
    user). The latter error value refers to error that are so serious
    that calculations can not be performed reliably (usually invalid
    input-files or arguments), or for some of the getter-functions,
    the calculation hasn't been performed yet. The return value
    SASALIB_SUCCESS means no errors have been spotted.

    Functions that have real-valued return values return negative
    numbers if calculation failed for some reason. The documentation
    for each function explains when this can happen.

*/

#include <stdio.h>

//#include "structure.h"

typedef struct sasalib_ sasalib_t;
typedef enum {SASALIB_LEE_RICHARDS, SASALIB_SHRAKE_RUPLEY} 
    sasalib_algorithm;

/** 4 classes of atoms/chemical groups used */
typedef enum {
    SASALIB_POLAR=0, SASALIB_APOLAR, 
    SASALIB_NUCLEICACID, SASALIB_CLASS_UNKNOWN
} sasalib_class;

// Limit for protein name lengths
#define SASALIB_NAME_LIMIT 30

// Default parameters
#define SASALIB_DEF_PROBE_RADIUS 1.4
#define SASALIB_DEF_SR_N 100
#define SASALIB_DEF_LR_D 0.25

#define SASALIB_SUCCESS 0
#define SASALIB_FAIL -1 /* Errors that mean results will be undefined */
#define SASALIB_WARN -2 /* Errors that mean calculation can proceed
			 * but will not comply to specification. */

/*****************************************/
/** Initialization, freeing and copying **/
/*****************************************/

/** Allocates empty sasalib_t object with default parameters. Must be
    called before calculation. */
sasalib_t* sasalib_init();

/** Destroys s. If s was initialized directly from a PDB-file, the
    internally stored protein will be freed. If the protein was
    supplied by the user, it will not be freed by this function. */
void sasalib_free(sasalib_t *s);

/** Copies all algorithm parameters of source to target. Data from
    calculations are not copied, i.e. this can be used for setting up
    several calculations with identical parameters for different
    proteins. */
void sasalib_copy_param(sasalib_t *target, const sasalib_t *source);

/**************************/
/** Perform calculations **/
/**************************/

/** Calculate SASA from given atom coordinates and radii. Doesn't even
    have to be a protein. */
int sasalib_calc_coord(sasalib_t*, const double *coord, 
                       const double *r, size_t n);

/** performs calculations on PDB-file. Results stored in parameter
    s. If s is not initialized default values are used, these are
    stored in s. Returns SASALIB_SUCCESS if calculation successful,
    prints an error and returns SASALIB_FAIL if not. If the object s
    has been used in calculations previously, the results from these
    will be over-written. */
int sasalib_calc_pdb(sasalib_t *s, FILE *pdb_file);

/** Reads pdb-file and calculates radii for each atom. Memory is
    allocated to store them in the array 'r'. The return value is the
    size of the allocated array. Prints error and returns SASALIB_FAIL
    if reading input fails. This can be used if coordinates are linked
    as below and default radii are to be used. 
    Not properly tested yet!
*/
int sasalib_generate_radii(double **r, FILE *pdb_file);

/** Link a set of coordinates to the sasalib_t object, if these
    coordinates are updated sasalib_refresh(1) can be used to
    recalculate SASA. Sasalib will not change the coordinates.  The
    array coord is of size 3*n with coordinates in the order
    x1,y1,z1,x2,y2,z2,... The array r should containt the radius of
    each atom and be of size n. If the sasalib_t-object has been
    initalized with a PDB file or structure_t-object these are
    released. */
int sasalib_link_coord(sasalib_t*, const double *coord,
                       double *r, size_t n);

/** Recalculates SASA, based on the assumption that a set of external
    coordinates have been updated elsewhere. Returns SASALIB_FAIL if
    no coordinates or radii are found in s. Returns SASALIB_SUCCESS
    upon successful computation. */
int sasalib_refresh(sasalib_t *s);

/** Returns the number of atoms in the latest SASA
    calculation. Returns 0 if no coordinates have been linked or no
    calculation has been performed. */
size_t sasalib_n_atoms(const sasalib_t*);

/******************************/
/** Settings for calculation **/
/******************************/

/** Sets algorithm. Returns SASALIB_SUCCESS if alg is valid, returns
    SASALIB_WARN else. */
int sasalib_set_algorithm(sasalib_t*, sasalib_algorithm alg);

/** Returns algorithm. */
sasalib_algorithm sasalib_get_algorithm(const sasalib_t*);

/** Returns name of algorithm. */
const char* sasalib_algorithm_name(const sasalib_t*);

/** Sets probe-radius for SASA calculations (default 1.4 Å). Returns
    SASALIB_SUCCESS for valid r-values. Prints error message, returns
    SASALIB_WARN and uses default (SASALIB_DEF_PROBE_RADIUS) else. */
int sasalib_set_probe_radius(sasalib_t*,double r);

/** Returns probe radius. */
double sasalib_get_probe_radius(const sasalib_t*);

/** Sets number of points for S&R algorithm (default 100). Returns
    SASALIB_SUCCESS if n is valid, prints error message, returns
    SASALIB_WARN and sets to default value (SASALIB_DEF_SR_N) else. */
int sasalib_set_sr_points(sasalib_t*, int n);

/** Returns number of points for S&R algorithm. Returns SASALIB_WARN
    if S&R algorithm not selected. */
int sasalib_get_sr_points(const sasalib_t*);

/** Sets slice width d for L&R algorithm in Ångström (default 0.25
    Å). Returns SASALIB_SUCCESS if d is valid, prints error message,
    returns SASALIB_WARN and sets to default value
    (SASALIB_DEF_LR_D) else. */
int sasalib_set_lr_delta(sasalib_t*, double d);

/** Returns slice width for L&R algorithm in Ångström. Returns negative
    value if L&R algorithm not selected. */
double sasalib_get_lr_delta(const sasalib_t*);

/** Sets the number of threads for parallel computation, useful for
    large proteins and high resolution. Returns SASALIB_SUCCESS if n
    is valid, uses default value and returns SASALIB_WARN else. */
int sasalib_set_nthreads(sasalib_t*,int n);

/** Returns the number of threads used in the calcultion */
int sasalib_get_nthreads(const sasalib_t*);

/** Sets name of protein, useful for logging. Uses last
    SASALIB_NAME_LIMIT characters if name is too long (since then it
    is probably a file name and the last characters are the most
    interesting). */
void sasalib_set_proteinname(sasalib_t*,const char*);

/** Returns protein name. */
const char* sasalib_get_proteinname(const sasalib_t*);

/*************/
/** Results **/
/*************/

/** Returns the total SASA. Negative return value and warning printed
    if calculation hasn't been performed yet. */
double sasalib_area_total(const sasalib_t*);

/** Returns the SASA of class c
    (polar/apolar/nucleic/unknown). Returns negative value if
    calculation has not been performed yet. */
double sasalib_area_class(const sasalib_t*, sasalib_class c);

/** Prints name/value-pairs with the total SASA of each residue
    type. The standard 20 amino acids are always included in output,
    non-standard ones and nucleotides only if they were present in
    input. Return SASALIB_FAIL if file-pointer NULL or if no
    calculation has been performed yet.*/
int sasalib_per_residue(FILE *output, const sasalib_t *s);

/** Returns total SASA for residue of specified type. If the value of
    res_name is unknown, a warning is printed and the SASA value for
    residue type UNK returned, because that is where it would be
    stored. I.e. if residues not known by Sasalib are used, the only
    option currently is to group them under "UNK". Returns negative
    value if called before calculations have been performed.*/
double sasalib_area_residue(const sasalib_t*, const char *res_name);

/** Takes original PDB and replaces B-factors with those from latest
    calculation. Returns SASALIB_FAIL if there is no previous PDB
    input to base output on, if there are problems with the output
    destination, if there are no SASA-values to use, or there are
    inconsistencies between stored structure and SASA-values. */
int sasalib_write_pdb(FILE *output, const sasalib_t *s);

/**********************************/
/** Results for individual atoms **/
/**********************************/

/** Returns SASA value for atom i. Prints error and returns negative
    value if atom index is invalid or if no calculation has been
    performed. */
double sasalib_area_atom(const sasalib_t*, int i);

/** Returns array containg SASA for all atoms. Returns NULL if no
    results available. */
const double* sasalib_area_atom_array(const sasalib_t*);

/** Returns radius of atom i. Prints error and returns negative
    value if atom index is invalid. */
double sasalib_radius_atom(const sasalib_t*, int i);

/** Returns array containing all atomic radii. Returns NULL if no results
    are available. */
const double* sasalib_radius_atom_array(const sasalib_t *s);
/****************************/
/** Other types of results **/
/****************************/

/** Prints log of calculation results to specified file. Returns
    SASALIB_SUCCESS on success, SASALIB_WARN if inconsistencies are
    detected (with explanatory error-message). */
int sasalib_log(FILE *log, const sasalib_t*);


#endif
