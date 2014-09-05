/*
  Copyright Simon Mitternacht 2013-2014.

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

#ifndef FREESASA_H
#define FREESASA_H

/** freesasa.h provides functions to init and perform a SASA
    calculation using FreeSASA with standard atom classifications. The
    user may optionally select algorithm and provide parameters.

    * Coordinates *

    If users wish to supply their own coordinates and radii, these are
    accepted as arrays of doubles. The coordinate-array should have
    size 3*n with coordinates in the order x1,y1,z1,x2,y2,z2,...

    * Error-reporting *

    All errors are written to stderr and are prefixed with the string
    'freesasa'. There are two error return values FREESASA_WARN and
    FREESASA_FAIL (see documentation of each function to see when these
    are used). The former refers to minor errors where a calculations
    can still proceed, but details might not be as intended
    (e.g. default parameters are used instead of those provided by
    user). The latter error value refers to error that are so serious
    that calculations can not be performed reliably (usually invalid
    input-files or arguments), or for some of the getter-functions,
    the calculation hasn't been performed yet. The return value
    FREESASA_SUCCESS means no errors have been spotted.

    Functions that have real-valued return values return negative
    numbers if calculation failed for some reason. The documentation
    for each function explains when this can happen.

*/

#include <stdio.h>

//#include "structure.h"

typedef struct freesasa_ freesasa_t;
typedef enum {FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY}
    freesasa_algorithm;

/** 4 classes of atoms/chemical groups used */
typedef enum {
    FREESASA_POLAR=0, FREESASA_APOLAR,
    FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN
} freesasa_class;

// Limit for protein name lengths
#define FREESASA_NAME_LIMIT 30

// Default parameters
#define FREESASA_DEF_PROBE_RADIUS 1.4
#define FREESASA_DEF_SR_N 100
#define FREESASA_DEF_LR_D 0.25

#define FREESASA_SUCCESS 0
#define FREESASA_FAIL -1 /* Errors that mean results will be undefined */
#define FREESASA_WARN -2 /* Errors that mean calculation can proceed
                          * but will not comply to specification. */

/*****************************************/
/** Initialization, freeing and copying **/
/*****************************************/

#ifdef __cplusplus
extern "C"{
#endif

/** Allocates empty freesasa_t object with default parameters. Must be
    called before calculation. */
    freesasa_t* freesasa_init();

/** Destroys s. If s was initialized directly from a PDB-file, the
    internally stored protein will be freed. If the protein was
    supplied by the user, it will not be freed by this function. */
    void freesasa_free(freesasa_t *s);

/** Copies all algorithm parameters of source to target. Data from
    calculations are not copied, i.e. this can be used for setting up
    several calculations with identical parameters for different
    proteins. */
    void freesasa_copy_param(freesasa_t *target, const freesasa_t *source);

/**************************/
/** Perform calculations **/
/**************************/

/** Calculate SASA from given atom coordinates and radii. Doesn't even
    have to be a protein. Coordinates and radii will not be
    stored. Return FREESASA_SUCESS upon successful calculation,
    prints and error and returns FREESASA_FAIL else. */
    int freesasa_calc_coord(freesasa_t*, const double *coord,
                            const double *r, size_t n);

/** performs calculations on PDB-file. Results stored in parameter
    s. If s is not initialized default values are used, these are
    stored in s. Returns FREESASA_SUCCESS if calculation successful,
    prints an error and returns FREESASA_FAIL if not. If the object
    has been used in calculations previously, the results from these
    will be over-written. */
    int freesasa_calc_pdb(freesasa_t *s, FILE *pdb_file);

/** Reads pdb-file and calculates radii for each atom. Memory is
    allocated to store them in the array 'r'. The return value is the
    size of the allocated array. Prints error and returns FREESASA_FAIL
    if reading input fails. This can be used if coordinates are linked
    as below and default radii are to be used.
    Not properly tested yet!
*/
//int freesasa_generate_radii(double **r, FILE *pdb_file);

/** Returns default radius of an atom type. The residue and atom names
    are the default names used in PDB files. Any whitespace in the
    standard needs to be included here, i.e. "CA" should be " CA ".
    Unknown atom types and hydrogens are assigned radius 0.0. */
    double freesasa_radius(const char* residue_name, const char* atom_name);

/** Link a set of coordinates to the freesasa_t object, if these
    coordinates are updated freesasa_refresh(1) can be used to
    recalculate SASA. FreeSASA will not change the coordinates.  The
    array coord is of size 3*n with coordinates in the order
    x1,y1,z1,x2,y2,z2,... The array r should containt the radius of
    each atom and be of size n. If the freesasa_t-object has been
    initalized with a PDB file or structure_t-object these are
    released. */
    int freesasa_link_coord(freesasa_t*, const double *coord,
                            double *r, size_t n);

/** Recalculates SASA, based on the assumption that a set of external
    coordinates have been updated elsewhere. Returns FREESASA_FAIL if
    no coordinates or radii are found in s. Returns FREESASA_SUCCESS
    upon successful computation. */
    int freesasa_refresh(freesasa_t *s);

/** Returns the number of atoms in the latest SASA
    calculation. Returns 0 if no coordinates have been linked or no
    calculation has been performed. */
    size_t freesasa_n_atoms(const freesasa_t*);

/******************************/
/** Settings for calculation **/
/******************************/

/** Sets algorithm. Returns FREESASA_SUCCESS if alg is valid, returns
    FREESASA_WARN else. */
    int freesasa_set_algorithm(freesasa_t*, freesasa_algorithm alg);

/** Returns algorithm. */
    freesasa_algorithm freesasa_get_algorithm(const freesasa_t*);

/** Returns name of algorithm. */
    const char* freesasa_algorithm_name(const freesasa_t*);

/** Sets probe-radius for SASA calculations (default 1.4 Å). Returns
    FREESASA_SUCCESS for valid r-values. Prints error message, returns
    FREESASA_WARN and uses default (FREESASA_DEF_PROBE_RADIUS) else. */
    int freesasa_set_probe_radius(freesasa_t*,double r);

/** Returns probe radius. */
    double freesasa_get_probe_radius(const freesasa_t*);

/** Sets number of points for S&R algorithm (default 100). Returns
    FREESASA_SUCCESS if n is valid, prints error message, returns
    FREESASA_WARN and sets to default value (FREESASA_DEF_SR_N) else. */
    int freesasa_set_sr_points(freesasa_t*, int n);

/** Returns number of points for S&R algorithm. Returns FREESASA_WARN
    if S&R algorithm not selected. */
    int freesasa_get_sr_points(const freesasa_t*);

/** Sets slice width d for L&R algorithm in Ångström (default 0.25
    Å). Returns FREESASA_SUCCESS if d is valid, prints error message,
    returns FREESASA_WARN and sets to default value
    (FREESASA_DEF_LR_D) else. */
    int freesasa_set_lr_delta(freesasa_t*, double d);

/** Returns slice width for L&R algorithm in Ångström. Returns negative
    value if L&R algorithm not selected. */
    double freesasa_get_lr_delta(const freesasa_t*);

/** Sets the number of threads for parallel computation, useful for
    large proteins and high resolution. Returns FREESASA_SUCCESS if n
    is valid, uses default value and returns FREESASA_WARN else. */
    int freesasa_set_nthreads(freesasa_t*,int n);

/** Returns the number of threads used in the calcultion */
    int freesasa_get_nthreads(const freesasa_t*);

/** Sets name of protein, useful for logging. Uses last
    FREESASA_NAME_LIMIT characters if name is too long (since then it
    is probably a file name and the last characters are the most
    interesting). */
    void freesasa_set_proteinname(freesasa_t*,const char*);

/** Returns protein name. */
    const char* freesasa_get_proteinname(const freesasa_t*);

/*************/
/** Results **/
/*************/

/** Returns the total SASA. Negative return value and warning printed
    if calculation hasn't been performed yet. */
    double freesasa_area_total(const freesasa_t*);

/** Returns the SASA of class c
    (polar/apolar/nucleic/unknown). Returns negative value if
    calculation has not been performed yet. */
    double freesasa_area_class(const freesasa_t*, freesasa_class c);

/** Prints name/value-pairs with the total SASA of each residue
    type. The standard 20 amino acids are always included in output,
    non-standard ones and nucleotides only if they were present in
    input. Return FREESASA_FAIL if file-pointer NULL or if no
    calculation has been performed yet.*/
    int freesasa_per_residue(FILE *output, const freesasa_t *s);

/** Returns total SASA for residue of specified type. If the value of
    res_name is unknown, a warning is printed and the SASA value for
    residue type UNK returned, because that is where it would be
    stored. I.e. if residues not known by FreeSASA are used, the only
    option currently is to group them under "UNK". Returns negative
    value if called before calculations have been performed.*/
    double freesasa_area_residue(const freesasa_t*, const char *res_name);

/** Takes original PDB and replaces B-factors with those from latest
    calculation. Returns FREESASA_FAIL if there is no previous PDB
    input to base output on, if there are problems with the output
    destination, if there are no SASA-values to use, or there are
    inconsistencies between stored structure and SASA-values. */
    int freesasa_write_pdb(FILE *output, const freesasa_t *s);

/**********************************/
/** Results for individual atoms **/
/**********************************/

/** Returns SASA value for atom i. Prints error and returns negative
    value if atom index is invalid or if no calculation has been
    performed. */
    double freesasa_area_atom(const freesasa_t*, int i);

/** Returns array containg SASA for all atoms. Returns NULL if no
    results available. */
    const double* freesasa_area_atom_array(const freesasa_t*);

/** Returns radius of atom i. Prints error and returns negative
    value if atom index is invalid or no value available. */
    double freesasa_radius_atom(const freesasa_t*, int i);

/** Returns array containing all atomic radii. Returns NULL if no results
    are available. */
    const double* freesasa_radius_atom_array(const freesasa_t *s);
/****************************/
/** Other types of results **/
/****************************/

/** Prints log of calculation results to specified file. Returns
    FREESASA_SUCCESS on success, FREESASA_WARN if inconsistencies are
    detected (with explanatory error-message). */
    int freesasa_log(FILE *log, const freesasa_t*);

#ifdef __cplusplus
}
#endif

#endif
