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

#ifndef SASALIB_H
#define SASALIB_H

/** sasalib.h provides functions to init and perform a SASA
    calculation using Sasalib with standard atom classifications. The
    user may optionally select algorithm and provide parameters. */

#include <stdio.h>
#include "structure.h"

typedef struct sasalib_ sasalib_t;
typedef enum {LEE_RICHARDS, SHRAKE_RUPLEY} sasalib_algorithm;

// Limit for protein name lengths
#define SASALIB_NAME_LIMIT 30

// Default parameters
#define SASALIB_DEF_SR_N 100
#define SASALIB_DEF_LR_D 0.25


/** Allocates empty sasalib_t obejct with default parameters. Must be
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

/** performs calculations on PDB-file. Results stored in parameter
    s. If s is not initialized default values are used, these are
    stored in s. Returns 0 if calculation successful, prints an error
    and returns 1 if not.*/
int sasalib_calc_pdb(sasalib_t *s, FILE *pdb_file);

/** Same as sasalib_calcpdb but from already initialized protein. */
int sasalib_calc_protein(sasalib_t*, structure_t*);

/** Sets algorithm. Returns 0 if alg is valid, returns 1 else. */
int sasalib_set_algorithm(sasalib_t*, sasalib_algorithm alg);

/** Returns algorithm. */
sasalib_algorithm sasalib_get_algorithm(const sasalib_t*);

/** Returns name of algorithm. */
const char* sasalib_algorithm_name(const sasalib_t*);

/** Sets number of points for S&R algorithm (default 100). Returns 0
    if n is valid, prints error message and returns 1 else. */
int sasalib_set_sr_points(sasalib_t*, int n);

/** Returns number of points for S&R algorithm. Returns negative
    number if S&R algorithm not selected. */
int sasalib_get_sr_points(const sasalib_t*);

/** Sets slice width for L&R algorithm in Ångström (default 0.25
    Å). Returns 0 if n is valid, prints error message and returns 1
    else. */
int sasalib_set_lr_delta(sasalib_t*, double d);

/** Returns slice width for L&R algorithm in Ånström. Returns negative
    value if L&R algorithm not selected. */
int sasalib_get_lr_delta(const sasalib_t*);

#ifdef PTHREADS
/** Sets the number of threads for parallel computation, useful for
    large proteins and high resolution. Returns 0 if n is valid, 1
    else. */
int sasalib_set_nthreads(sasalib_t*,int n);

/** Returns the number of threads used in the calcultion */
int sasalib_get_nthreads(const sasalib_t*);
#endif

/** Sets name of protein, useful for logging. Uses last
    SASALIB_NAME_LIMIT characters if name is too long (since then it
    is probably a file name and the last characters are the most
    interesting). */
void sasalib_set_proteinname(sasalib_t*,const char*);

/** Returns protein name. */
const char* sasalib_get_proteinname(const sasalib_t*);

/** Returns the total SASA. Negative return value and warning printed
    if calculation hasn't been performed yet. */
// not implemented yet
double sasalib_area_total(const sasalib_t*);

/** Returns the polar SASA. Negative return value and warning printed
    if calculation hasn't been performed yet. */
double sasalib_area_polar(const sasalib_t*);

/** Returns the apolar SASA. Negative return value and warning printed
    if calculation hasn't been performed yet. */
double sasalib_area_apolar(const sasalib_t*);

/** Returns the SASA of any nucleic acids present in the
    structure. Negative return value and warning printed if
    calculation hasn't been performed yet. */
double sasalib_area_nucleicacid(const sasalib_t*);

/** Returns SASA value for atom i. Prints error and returns negative
    value if atom index is invalid. */
double sasalib_area_atom(const sasalib_t*, int i);

/** Returns radius of atom i. Prints error and returns negative
    value if atom index is invalid. */
double sasalib_radius_atom(const sasalib_t*, int i);

/** Prints log to specified file (after the fact). Returns 0 on
    success, 1 on failure. */
int sasalib_log(FILE *log, const sasalib_t*);

/** Prints the total SASA for each residue of the structure */
void sasalib_per_residue(FILE *output, const sasalib_t*);

#endif
