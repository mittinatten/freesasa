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

#ifndef FREESASA_H
#define FREESASA_H

/**
    @file
    @author Simon Mitternacht
    
    @section The FreeSASA API

    The header freesasa.h contains the API of FreeSASA and is the only
    header installed by the `make install` target. It provides
    functions to init and perform a SASA calculation. The user xan
    select algorithm and provide parameters. The type ::freesasa is
    used to store parameters and access results. The manual in
    `doc/manual.pdf` gives some examples of how to use this API.

    @subsection Coordinates

    If users wish to supply their own coordinates and radii, these are
    accepted as arrays of doubles. The coordinate-array should have
    size 3*n with coordinates in the order x1,y1,z1,x2,y2,z2,...

    @subsection Error-reporting 

    Input and user parameters are checked for errors and
    inconsistencies. All errors are written to stderr and are prefixed
    with the string 'freesasa'. There are two error codes
    ::FREESASA_WARN and ::FREESASA_FAIL (see documentation of each
    function to see when these are used). ::FREESASA_SUCCESS is used
    for success.

    Errors that are attributable to programmers using the library,
    such as passing null pointers, or calling functions in the wrong
    order, are checked by asserts.

    Memory allocation errors are only checked with asserts. These
    should be rare in a library of this type, the asserts are there to
    allow debugging should they occur.

    @subsection Thread-safety 
    
    The only state the library stores is the verbosity level (set by
    freesasa_set_verbosity()). It should be clear from the
    documentation when the other functions have side effects such as
    memory allocation and I/O, and thread-safety should generally not
    be an issue (to the extent that your c library has a threadsafe
    fprintf). The SASA calculation itself can be parallelized by
    increasing the number of threads through freesasa_set_nthreads()
    before calling any of the calculation-functions.
 */

#include <stdio.h>

//#ifdef __cplusplus
//extern "C"{
//#endif

/**
    The FreeSASA algorithms.
 */
typedef enum {FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY}
    freesasa_algorithm;

/**
    4 classes of atoms/chemical groups used 
 */
typedef enum {
    FREESASA_POLAR=0, FREESASA_APOLAR,
    FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN
} freesasa_class;

/**
   Granularity levels for arrays of results in freesasa_string_value_pairs().
 */
typedef enum  {FREESASA_ATOMS, FREESASA_RESIDUES}
    freesasa_result_type;

/**
    Verbosity levels. 
    - FREESASA_V_NORMAL: print all errors and warnings.
    - FREESASA_V_NOWARNINGS: print only errors.
    - FREESASA_V_SILENT: print no errors and warnings.
 */
typedef enum {FREESASA_V_NORMAL,
              FREESASA_V_NOWARNINGS,
              FREESASA_V_SILENT} freesasa_verbosity;

/// Limit for protein name lengths
#define FREESASA_NAME_LIMIT 30

// Default parameters
#define FREESASA_DEF_PROBE_RADIUS 1.4 //!< Default probe radius (in Ångtström)
#define FREESASA_DEF_SR_N 100 //!< Default number of test points in S&R
#define FREESASA_DEF_LR_D 0.25 //!< Default slice width in L&R (in Ångström)

#define FREESASA_SUCCESS 0 //!< All is ok
#define FREESASA_FAIL -1 //!< Something went seriously wrong.
#define FREESASA_WARN -2 //!< Something went wrong, but results might still be meaningful

//! Struct to store parameters for SASA calculation
typedef struct {
    freesasa_algorithm alg;
    double probe_radius;
    int shrake_rupley_n_points;
    double lee_richards_delta;
    int n_threads;
    int include_hetatm;
} freesasa_parameters;

//! Struct to store results of SASA calculation
typedef struct {
    double total;
    double *sasa;
    int n_atoms;
} freesasa_result;

typedef struct freesasa_structure freesasa_structure;

/**
    Struct used to store n string-value-pairs (strvp) in arrays of
    doubles and strings. freesasa_strvp_free() assumes both arrays
    and strings are dynamically allocated.
 */
typedef struct {
    double *value;
    char **string;
    int n;
} freesasa_strvp;

/**
    Type used to store SASA values for a groups of atoms and strings
    describing those groups of atoms.
 */
typedef freesasa_strvp freesasa_per_class;

struct freesasa_classifier;
/**
    Struct used for calculating classes and radii
 */
typedef struct freesasa_classifier {
    double (*radius)(const char* res_name,
                     const char* atom_name,
                     const struct freesasa_classifier*);
    int (*sasa_class)(const char* res_name,
                      const char* atom_name,
                      const struct freesasa_classifier*);
    const char* (*class2str)(int the_class,
                             const struct freesasa_classifier*);
    int n_classes;
    void *classifier_data;
} freesasa_classifier;

extern const freesasa_parameters freesasa_default_parameters;

extern const freesasa_classifier freesasa_default_classifier;

extern const freesasa_classifier freesasa_residue_classifier;

/**
    Calculates SASA based for a given structure and atomic radii.

    @param result Results are written to this pointer
    @param structure The structure
    @param radii Atomic radii, this array should have same number of
    elements as there are atoms in the structure.
    @param parameters Parameters for the calculation, if NULL
    defaults are used.

    @return ::FREESASA_SUCCESS if calculations were successful. Else
    ::FREESASA_FAIL
*/
int freesasa_calc_structure(freesasa_result *result,
                            const freesasa_structure *structure,
                            const double *radii,
                            const freesasa_parameters *parameters);

/**
    Calculates SASA based for a given set of coordinates and radii.

    @param result Results are written to this pointer.
    @param parameters Parameters for the calculation, if NULL
    defaults are used.
    @param xyz Array of coordinates in the form x1,y1,z1,x2,y2,z2,...,xn,yn,zn.
    @param radii Radii, this array should have same number of
    elements as there are coordinates.
    @param n Number of coordinates (i.e. xyz has size 3*n, radii size n).

    @return ::FREESASA_SUCCESS if calculations were successful. Else
    ::FREESASA_FAIL
 */
int freesasa_calc_coord(freesasa_result *result,
                        const double *xyz, 
                        const double *radii,
                        int n,
                        const freesasa_parameters *parameters);

/**
    Frees the contents of ::freesasa_result object.

    @param result the object to be freed.
 */
void freesasa_result_free(freesasa_result result);

freesasa_structure* freesasa_structure_from_pdb(FILE *pdb,
                                                int include_hetatm);

void freesasa_structure_free(freesasa_structure* structure);

double* freesasa_structure_radius(freesasa_structure *structure,
                                  freesasa_classifier *classifier);
/**
    Generate a classifier from a config-file.

    @param file File containing configuration
    @return The generated classifier. NULL if file there were
    problems parsing or reading the file.
 */
freesasa_classifier* freesasa_classifier_from_file(FILE *file);

/**
    Frees the contents of a classifier object

    @param classifier The classifier.
 */
void freesasa_classifier_free(freesasa_classifier classifier);

/**
    Sums up the SASA for groups of atoms defined by a classifier.

    @param result The results to be analyzed.
    @param structure Structure to be used to determine atom types.
    @param classifier The classifier. If NULL, default is used.
    @return A new set of string-value-pairs if classifications was
    successful that should be freed with
    freesasa_strvp_free(). Returns NULL if classifier was not
    compatible with structure.
 */
freesasa_strvp* freesasa_result_classify(freesasa_result result,
                                         const freesasa_structure *structure,
                                         const freesasa_classifier *classifier);

/**
    Frees a ::freesasa_strvp object

    @param strvp the object to be freed
 */
void freesasa_strvp_free(freesasa_strvp *strvp);

/**
    Write SASA values and atomic radii to PDB-file.

    Takes original PDB and for each atom replace the B-factor with the
    atom's SASA values from latest calculation and the occupancy
    factor with the radius used.

    @param output File to write to.
    @param result SASA values.
    @param structure Structure to use to print PDB.
    @param radii Radii of atoms.
    @return ::FREESASA_FAIL if there is no previous PDB input to base
    output on. ::FREESASA_SUCCESS else.
 */
int freesasa_write_pdb(FILE *output, 
                       freesasa_result result,
                       const freesasa_structure *structure, 
                       const double *radii);

/**
    Print SASA for all residue types to file.

    Prints name/value-pairs with the total SASA of each residue
    type. The standard 20 amino acids are always included in output,
    non-standard ones and nucleotides only if they were present in
    input. Each line in the output is prefixed by the string 'RES:'.

    @param output Output file.
    @param result SASA values.
    @param structure The structure (includes sequence information).
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.
*/
int freesasa_per_residue_type(FILE *output, 
                              freesasa_result result,
                              const freesasa_structure *structure);
/**
    Print SASA for each residue individually to file. 

    Each line in the output is prefixed by the string 'SEQ:'.

    @param output Output file.
    @param result SASA values.
    @param structure The structure (includes sequence information).
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.
*/
int freesasa_per_residue(FILE *output,
                         freesasa_result result,
                         const freesasa_structure *structure);
/**
    Log calculation results.

    Prints log of calculation results to specified file. 

    @param log Output-file.
    @param result SASA values.
    @param name Name of the protein
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if
    problems writing to file.
*/
int freesasa_log(FILE *log, 
                 freesasa_result result,
                 const char *name,
                 const freesasa_structure *structure, 
                 const freesasa_parameters *parameters,
                 const freesasa_strvp* area_of_classes);

/**
    Set the global verbosity level.

    @arg v the verbosity level
    @return ::FREESASA_SUCCESS. If v is invalid ::FREESASA_FAIL.
    @see freesasa_verbosity
*/
int freesasa_set_verbosity(freesasa_verbosity v);

/**
    Get the current verbosity level

    @return the verbosity level. 
*/
freesasa_verbosity freesasa_get_verbosity(void);

//#ifdef __cplusplus
//}
//#endif

#endif
