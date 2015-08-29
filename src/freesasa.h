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
    @defgroup API Public API
 */

/**
    @defgroup StructureAPI Structure API
    @ingroup API
 */

/**
    @file
    @author Simon Mitternacht
    @ingroup API

    Contains the @ref API of FreeSASA.
 */

/**
    @addtogroup API
    @section The FreeSASA API

    The header @ref freesasa.h contains the @ref API of FreeSASA and
    is the only header installed by the `make install` target. It
    provides the functions and data types necessary to perfrom and
    analyze a SASA calculation, including facilities to customize
    assignment of radii to, and classification of atoms. The API for
    dealing with structures is documented in the submodule @ref
    StructureAPI (but the functions are in the same header).

    @subsection Customizing Customizing behavior

    The types ::freesasa_parameters and ::freesasa_classifier can be
    used to change the parameters of the calculations. Users who wish
    to use the defaults can pass NULL wherever pointers to these are
    requested.

    @subsubsection Parameters Parameters

    Changing parameters is done by passing a ::freesasa_parameters
    object with the desired values. It can be initialized to default
    by

    ~~~{.c}
    freesasa_parameters p = freesasa_default_parameters;
    ~~~

    To allow the user to only change the parameters that are
    non-default.

    @subsubsection Classification Specifying atomic radii and classes

    The type ::freesasa_classifier has function pointers to functions
    that take residue and atom names as argument (pairs such as
    "ALA"," CA "), and returns a radius or a class (polar, apolar,
    etc). Such a classifier can be passed freesasa_structure_radius()
    to generate an array of atomic radii to be used to calculate the
    SASA of the structure. It can also be used in
    freesasa_result_classify() to get the SASA integrated over the
    different classes of atoms, i.e. the SASA of all polar atoms, etc.

    Users of the API can provide their own classification by writing
    their own functions and providing them via a ::freesasa_classifier
    object. A classifier-configuration can also be read from a file
    using freesasa_classifier_from_file() (see @ref Config-file).

    The default classifier is available as a global const variable
    ::freesasa_default_classifier.

    @subsubsection Config-file Classifier configuration files

    The configuration files read by freesasa_classifier_from_file()
    should have two sections: `types:` and `atoms:`.

    The types-section defines what types of atoms are available
    (aliphatic, aromatic, hydroxyl, ...), what the radius of that type
    is and what class a type belongs to (polar, apolar, ...). The
    types are just a shorthand to associate a certain atom with a
    certain class and radius. The user is free to define as many types
    and classes as necessary.

    The atoms-section consists of triplets of residue-name, atom-name
    (as in the corresponding PDB entries) and type. A prototype file
    would be

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

    @subsection Coordinates

    If users wish to supply their own coordinates and radii, these are
    accepted as arrays of doubles passed to the function
    freesasa_calc_coord(). The coordinate-array should have size 3*n
    with coordinates in the order `x1,y1,z1,x2,y2,z2,...,xn,yn,zn`.

    @subsection Error-reporting 

    Errors due to user or system errors, such as malformatted
    config-files, I/O errors are reported through return values,
    either ::FREESASA_FAIL or ::FREESASA_WARN, or by NULL
    pointers. See the documentation for the individual functions.

    Errors that are attributable to programmers using the library,
    such as passing null pointers are checked by asserts.

    Memory allocation errors are only checked with asserts. These
    should be rare in a library of this type, the asserts are there to
    allow debugging should they occur.

    @subsection Thread-safety 
    
    The only state the library stores is the verbosity level (set by
    freesasa_set_verbosity()). It should be clear from the
    documentation when the other functions have side effects such as
    memory allocation and I/O, and thread-safety should generally not
    be an issue (to the extent that your C library has a threadsafe
    fprintf). The SASA calculation itself can be parallelized by
    passing a ::freesasa_parameters struct with
    ::freesasa_parameters.n_threads set to a value > 1 to
    freesasa_calc_pdb() or freesasa_calc_coord().
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

//! The FreeSASA algorithms. @ingroup API
typedef enum {FREESASA_LEE_RICHARDS, FREESASA_SHRAKE_RUPLEY}
    freesasa_algorithm;

/**
    Verbosity levels. 
    - FREESASA_V_NORMAL: print all errors and warnings.
    - FREESASA_V_NOWARNINGS: print only errors.
    - FREESASA_V_SILENT: print no errors and warnings.

    @ingroup API
*/    
typedef enum {FREESASA_V_NORMAL,
              FREESASA_V_NOWARNINGS,
              FREESASA_V_SILENT} freesasa_verbosity;
/**
    4 classes of atoms/chemical groups 
    (classes in freesasa_default_classifier)
    @ingroup API
 */
typedef enum {
    FREESASA_POLAR=0, FREESASA_APOLAR,
    FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN
} freesasa_class;


// Default parameters
#define FREESASA_DEF_PROBE_RADIUS 1.4 //!< Default probe radius (in Ångström) @ingroup API
#define FREESASA_DEF_SR_N 100 //!< Default number of test points in S&R @ingroup API
#define FREESASA_DEF_LR_D 0.25 //!< Default slice width in L&R (in Ångström) @ingroup API

#define FREESASA_SUCCESS 0 //!< All is ok @ingroup API
#define FREESASA_FAIL -1 //!< Something went seriously wrong. @ingroup API
#define FREESASA_WARN -2 //!< Something went wrong, but results might still be meaningful @ingroup API

//! Struct to store parameters for SASA calculation @ingroup API
typedef struct {
    freesasa_algorithm alg;     //!< Algorithm
    double probe_radius;        //!< Probe radius (in Ångström)
    int shrake_rupley_n_points; //!< Number of test points in S&R calculation
    double lee_richards_delta;  //!< Slice width in L&R calculation (in Ångström)
    int n_threads;              //!< Number of threads to use, if compiled with thread-support
} freesasa_parameters;

//! Struct to store results of SASA calculation @ingroup API
typedef struct {
    double total; //!< Total SASA in Ångström^2
    double *sasa; //!< SASA of each atom in Ångström^2
    int n_atoms;  //!< Number of atoms
} freesasa_result;

/**
    Struct for structure object.

    The struct includes coordinates, and atom names, etc. If it was
    initiated from a PDB file enough info will be stored so that
    a new PDB-file can be printed.

    @ingroup StructureAPI
*/
typedef struct freesasa_structure freesasa_structure;

/**
    Struct used to store n string-value-pairs (strvp) in arrays of
    doubles and strings. freesasa_strvp_free() assumes both arrays
    and strings are dynamically allocated.

    @ingroup API
 */
typedef struct {
    double *value; //!< Array of values
    char **string; //!< Array of strings
    int n;         //!< Number of values and strings
} freesasa_strvp;

/**
    Struct used for calculating classes and radii for atoms given
    their residue-names ('ALA','ARG',...) and atom-names
    ('CA','N',...).

    @ingroup API
 */
typedef struct freesasa_classifier {
    int n_classes; //!< Total number of different classes
    void *config;  //!< Optional configuration to allow flexibility

    //! Function that returns an atom radius.
    double (*radius)(const char* res_name,
                     const char* atom_name,
                     const struct freesasa_classifier*);

    //! Function that returns the class [0,1,...,n_classes-1] of an atom
    int (*sasa_class)(const char* res_name,
                      const char* atom_name,
                      const struct freesasa_classifier*);

    //! Function that converts a class to its string descriptor
    const char* (*class2str)(int the_class,
                             const struct freesasa_classifier*);

    //! Function that can be called to free the config-pointer
    void (*free_config)(void*);
} freesasa_classifier;

//! The default parameters for FreeSASA @ingroup API
extern const freesasa_parameters freesasa_default_parameters;

//! The default classifier, uses the functions in `classify.h` @ingroup API
extern const freesasa_classifier freesasa_default_classifier;
//! Classifier that classifies each atom according to residue @ingroup API
extern const freesasa_classifier freesasa_residue_classifier;

/**
    Calculates SASA based on a given structure and atomic radii.

    Return value is dynamically allocated, should be freed with
    freesasa_result_free().

    @param structure The structure
    @param radii Atomic radii, this array should have same number of
    elements as there are atoms in the structure.
    @param parameters Parameters for the calculation, if NULL
    defaults are used.

    @return The result of the calculation, NULL if something went wrong.

    @ingroup API
*/
freesasa_result* freesasa_calc_structure(const freesasa_structure *structure,
                                         const double *radii,
                                         const freesasa_parameters *parameters);

/**
    Calculates SASA based on a given set of coordinates and radii.

    Return value is dynamically allocated, should be freed with
    freesasa_result_free().

    @param parameters Parameters for the calculation, if NULL
    defaults are used.
    @param xyz Array of coordinates in the form x1,y1,z1,x2,y2,z2,...,xn,yn,zn.
    @param radii Radii, this array should have n elements..
    @param n Number of coordinates (i.e. xyz has size 3*n, radii size n).

    @return The result of the calculation, NULL if something went wrong.

    @ingroup API
 */
freesasa_result* freesasa_calc_coord(const double *xyz, 
                                     const double *radii,
                                     int n,
                                     const freesasa_parameters *parameters);

/**
    Frees a ::freesasa_result object.

    @param result the object to be freed.

    @ingroup API
 */
void freesasa_result_free(freesasa_result *result);

/**
    Init structure with coordinates from pdb-file.

    Reads in a PDB-file and generates a structure object.
    Automatically skips hydrogens. If an atom has alternative
    coordinates, only the first alternative is used. If a file has
    more than one `MODEL` (as in NMR structures) only the first model
    is used. User specifies if `HETATM` entries should be included. If
    non-default behavior is wanted, the PDB-file needs to be modified
    before calling this function, or atoms can be added manually one
    by one using freesasa_structure_add_atom().

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @param pdb Input PDB-file.
    @param include_hetatm The value 0 means only read `ATOM` entries, 1
    means also include `HETATM` entries.
    @return The generated structure. Returns `NULL` and prints error if
    input is invalid.

    @ingroup StructureAPI
*/
freesasa_structure* freesasa_structure_from_pdb(FILE *pdb,
                                                int include_hetatm);

/**
    Get number of atoms.
    
    @param s Self.
    @return Number of atoms.

    @ingroup StructureAPI
*/
int freesasa_structure_n(const freesasa_structure *s);

/**
    Free structure.

    @param structure The structure to free.

    @ingroup StructureAPI
 */
void freesasa_structure_free(freesasa_structure* structure);

/**
    Calculates radii to all atoms in the structure using provided
    classifier.

    Return value is dynamically allocated, should be freed with
    standard free().

    @param structure The structure.
    @param classifier The classifier. If NULL the default is used.
    @return Array of radii.

    @ingroup StructureAPI
 */
double* freesasa_structure_radius(const freesasa_structure *structure,
                                  const freesasa_classifier *classifier);
/**
    Generate a classifier from a config-file.

    File-format is documented at @ref Config-file

    Return value is dynamically allocated, should be freed with
    freesasa_classifier_free().

    @param file File containing configuration
    @return The generated classifier. NULL if file there were
    problems parsing or reading the file.

    @ingroup API
 */
freesasa_classifier* freesasa_classifier_from_file(FILE *file);

/**
    Frees the contents of a classifier object

    @param classifier The classifier.

    @ingroup API
 */
void freesasa_classifier_free(freesasa_classifier *classifier);

/**
    Sums up the SASA for groups of atoms defined by a classifier.

    Return value is dynamically allocated, should be freed with
    freesasa_strvp_free().

    @param result The results to be analyzed.
    @param structure Structure to be used to determine atom types.
    @param classifier The classifier. If NULL, default is used.
    @return A new set of string-value-pairs if classifications was
    successful. NULL if classifier was not compatible with structure.

    @ingroup API
 */
freesasa_strvp* freesasa_result_classify(const freesasa_result *result,
                                         const freesasa_structure *structure,
                                         const freesasa_classifier *classifier);

/**
    Frees a ::freesasa_strvp object

    @param strvp the object to be freed

    @ingroup API
 */
void freesasa_strvp_free(freesasa_strvp *strvp);

/**
    Write SASA values and atomic radii to new PDB-file.

    Takes PDB information from the provided structure and writes a new
    PDB-file to output, where the B-factors have been replaced with
    the atom's SASA values in the results, and the occupancy
    factors with the radii.

    @param output File to write to.
    @param result SASA values.
    @param structure Structure to use to print PDB.
    @param radii Radii of atoms.
    @return ::FREESASA_FAIL if there is no previous PDB input to base
    output on. ::FREESASA_SUCCESS else.

    @ingroup API
 */
int freesasa_write_pdb(FILE *output, 
                       freesasa_result *result,
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
    
    @ingroup API
*/
int freesasa_per_residue_type(FILE *output, 
                              freesasa_result *result,
                              const freesasa_structure *structure);

/**
    Print SASA for each residue individually to file. 

    Each line in the output is prefixed by the string 'SEQ:'.

    @param output Output file.
    @param result SASA values.
    @param structure The structure (includes sequence information).
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.

    @ingroup API
 */
int freesasa_per_residue(FILE *output,
                         freesasa_result *result,
                         const freesasa_structure *structure);

/**
    Log calculation results.

    Prints log of calculation results to specified file. 

    @param log Output-file.
    @param result SASA values.
    @param name Name of the protein, if NULL "unknown" used.
    @param parameters Parameters to print, if NULL defaults are used
    @param class_sasa The SASA values for each class, if NULL
    only total SASA printed
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if
    problems writing to file.

    @ingroup API
*/
int freesasa_log(FILE *log, 
                 freesasa_result *result,
                 const char *name,
                 const freesasa_parameters *parameters,
                 const freesasa_strvp* class_sasa);

/**
    Set the global verbosity level.

    @arg v the verbosity level
    @return ::FREESASA_SUCCESS. If v is invalid ::FREESASA_FAIL.
    @see freesasa_verbosity

    @ingroup API
*/
int freesasa_set_verbosity(freesasa_verbosity v);

/**
    Get the current verbosity level

    @return the verbosity level. 

    @ingroup API
 */
freesasa_verbosity freesasa_get_verbosity(void);

/**
    @addtogroup StructureAPI

    The Structure API contains functions to handle
    ::freesasa_structure objects. These can either be initialized from
    a PDB file (freesasa_structure_from_pdb()) or by adding atoms one
    by one using freesasa_structure_add_atom() (after initializing
    with freesasa_structure_new()). Structure objects are freed using
    freesasa_structure_free().

    The API also provides a number of functions to access properties
    of the structure and its atoms, which is not necessary for users
    of the core functionality of the library, but might be useful for
    extensions or more detailed analyses of the results.

    The Structure API is declared in the header @ref freesasa.h, just
    like the rest of the @ref API.
 */

/**
    Allocate empty structure.

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @return Empty structure.

    @ingroup StructureAPI
 */
freesasa_structure* freesasa_structure_new(void);

/**
    Add individual atom to structure.

    A structure can be built by adding atoms one by one. Storing
    residue numbers as strings allows for non-numeric labels. Will
    include hydrogens if added (i.e. up to caller to make sure these
    are excluded if necessesary).

    The three string arguments all have specific lengths, specified by
    the correspoding PDB ATOM fields. It is important to use the same
    padding as in the PDB specification (i.e. `" CA "` instead of
    `"CA"` or `" CA"`). (This might be more flexible in future
    versions of the library).

    @param s Self.
    @param atom_name String of 4 characters, of the format `" CA "`, `" OXT"`, etc.
    @param residue_name String of 3 charachters, of the format `"ALA"`, `"PHE"`, etc.
    @param residue_number String of 4 characters, of the format `"   1"`, `" 123"`, etc.
    @param chain_label Any character to label chain, typically `'A'`, `'B'`, etc.
    @param x x-coordinate of atom.
    @param y y-coordinate of atom.
    @param z z-coordinate of atom.
    @return ::FREESASA_SUCCESS if input valid. ::FREESASA_FAIL if any of
    the strings are malformatted. ::FREESASA_WARN if the atom type is
    unknown.

    @ingroup StructureAPI
 */
int freesasa_structure_add_atom(freesasa_structure *s,
                                const char* atom_name,
                                const char* residue_name,
                                const char* residue_number,
                                char chain_label,
                                double x, double y, double z);

/**
    Get atom name

    @param s Self.
    @param i Atom index.
    @return Atom name in the form `" CA "`, `" OXT"`, etc.

    @ingroup StructureAPI
 */
const char* freesasa_structure_atom_name(const freesasa_structure *s,
                                         int i);

/**
    Get residue name.

    @param s Self.
    @param i Atom index.
    @return Residue name in the form `"ALA"`, `"PHE"`, etc.

    @ingroup StructureAPI
 */
const char* freesasa_structure_atom_res_name(const freesasa_structure *s,
                                             int i);

/**
    Get residue number.

    @param s Self.
    @param i Atom index.
    @return Residue name in the form `"   1"`, `" 123"`, etc.

    @ingroup StructureAPI
 */
const char* freesasa_structure_atom_res_number(const freesasa_structure *s,
                                               int i);

/**
    Get chain label.

    @param s Self.
    @param i Atom index.
    @return Chain label (`'A'`, `'B'`, etc.)

    @ingroup StructureAPI
 */
char freesasa_structure_atom_chain(const freesasa_structure *s, int i);


#ifdef __cplusplus
}
#endif

#endif
