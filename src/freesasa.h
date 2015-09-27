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

    @brief Functions and datatypes for performing and analyzing SASA
    calculations.

    @section The FreeSASA API

    This header provides the functions and data types necessary to
    perfrom and analyze a SASA calculation using FreeSASA, including
    facilities to customize assignment of radii to, and classification
    of, atoms. There are also functions to access properties of a
    structure, to allow refined analysis of the results.
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
*/    
typedef enum {FREESASA_V_NORMAL,
              FREESASA_V_NOWARNINGS,
              FREESASA_V_SILENT} freesasa_verbosity;
/**
    4 classes of atoms/chemical groups 
    (classes in freesasa_default_classifier)
 */
typedef enum {
    FREESASA_POLAR=0, FREESASA_APOLAR,
    FREESASA_NUCLEICACID, FREESASA_CLASS_UNKNOWN
} freesasa_class;

// Default parameters
#define FREESASA_DEF_PROBE_RADIUS 1.4 //!< Default probe radius (in Ångström)
#define FREESASA_DEF_SR_N 100 //!< Default number of test points in S&R
#define FREESASA_DEF_LR_N 20 //!< Default number of slices per atom  in L&R
#ifdef HAVE_LIBPTHREAD
#define FREESASA_DEF_NUMBER_THREADS 2 //!< Default number of threads
#else
#define FREESASA_DEF_NUMBER_THREADS 1 //!< Default number of threads
#endif

// Important that success is 0 and failure is non-zero, don't change
#define FREESASA_SUCCESS 0 //!< All is ok
#define FREESASA_FAIL -1 //!< Something went seriously wrong. 
#define FREESASA_WARN -2 //!< Something went wrong, but results might still be meaningful.

// Parameters for reading structure from PDB
#define FREESASA_INCLUDE_HETATM 1 //!< Include HETATM entries
#define FREESASA_INCLUDE_HYDROGEN 2 //!< Include hydrogen atoms
#define FREESASA_SEPARATE_MODELS 4 //!< Read MODELs as separate structures
#define FREESASA_SEPARATE_CHAINS 8 //!< Read separate chains as separate structures
#define FREESASA_JOIN_MODELS 16 //!< Read MODELs as part of one big structure

//! Struct to store parameters for SASA calculation @ingroup API
typedef struct {
    freesasa_algorithm alg;       //!< Algorithm
    double probe_radius;          //!< Probe radius (in Ångström)
    int shrake_rupley_n_points;   //!< Number of test points in S&R calculation
    int lee_richards_n_slices;    //!< Number of slices per atom in L&R calculation
    int n_threads;                //!< Number of threads to use, if compiled with thread-support
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
    initiated from a PDB file enough info will be stored so that a
    PDB-file can be printed using the original one as template (see
    freesasa_write_pdb()).
 */
typedef struct freesasa_structure freesasa_structure;

/**
    Struct used to store n string-value-pairs (strvp) in arrays of
    doubles and strings. freesasa_strvp_free() assumes both arrays
    and strings are dynamically allocated.
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
 */
freesasa_result*
freesasa_calc_structure(const freesasa_structure *structure,
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
 */
freesasa_result*
freesasa_calc_coord(const double *xyz, 
                    const double *radii,
                    int n,
                    const freesasa_parameters *parameters);

/**
    Frees a ::freesasa_result object.

    @param result the object to be freed.
 */
void
freesasa_result_free(freesasa_result *result);

/**
    Generate a classifier from a config-file.

    Input file format described in @ref Config-file

    Return value is dynamically allocated, should be freed with
    freesasa_classifier_free().

    @param file File containing configuration 
    @return The generated classifier. NULL if there were problems
      parsing or reading the file or memory allocation problem.

    @see @ref Config-file
 */
freesasa_classifier*
freesasa_classifier_from_file(FILE *file);

/**
    Frees the contents of a classifier object

    @param classifier The classifier.
 */
void
freesasa_classifier_free(freesasa_classifier *classifier);

/**
    Sums up the SASA for groups of atoms defined by a classifier.

    Return value is dynamically allocated, should be freed with
    freesasa_strvp_free().

    @param result The results to be analyzed.
    @param structure Structure to be used to determine atom types.
    @param classifier The classifier. If NULL, default is used.
    @return A new set of string-value-pairs if classifications was
    successful. NULL if classifier was not compatible with structure
    or memory allocation failure.
 */
freesasa_strvp*
freesasa_result_classify(const freesasa_result *result,
                         const freesasa_structure *structure,
                         const freesasa_classifier *classifier);

/**
    Frees a ::freesasa_strvp object

    @param strvp the object to be freed
 */
void
freesasa_strvp_free(freesasa_strvp *strvp);

/**
    Write SASA values and atomic radii to new PDB-file.

    Takes PDB information from the provided structure and writes a new
    PDB-file to output, where the B-factors have been replaced with
    the atom's SASA values in the results, and the occupancy
    factors with the radii.

    Will only work if the structure was initialized from a PDB-file, i.e.
    using freesasa_structure_from_pdb().

    @param output File to write to.
    @param result SASA values.
    @param structure Structure to use to print PDB.
    @param radii Radii of atoms.

    @return ::FREESASA_SUCCESS if file written successfully.
    ::FREESASA_FAIL if there is no previous PDB input to base output
    on or if there were problems writing to the file.
 */
int
freesasa_write_pdb(FILE *output, 
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
 */
int
freesasa_per_residue_type(FILE *output, 
                          freesasa_result *result,
                          const freesasa_structure *structure);

/**
    Print SASA for each residue in the sequence to file.

    Each line in the output is prefixed by the string 'SEQ:'.

    @param output Output file.
    @param result SASA values.
    @param structure The structure (includes sequence information).
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.
 */
int
freesasa_per_residue(FILE *output,
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
 */
int
freesasa_log(FILE *log, 
             freesasa_result *result,
             const char *name,
             const freesasa_parameters *parameters,
             const freesasa_strvp* class_sasa);

/**
    Set the global verbosity level.

    @arg v the verbosity level
    @return ::FREESASA_SUCCESS. If v is invalid ::FREESASA_FAIL.
    @see freesasa_verbosity
 */
int
freesasa_set_verbosity(freesasa_verbosity v);

/**
    Get the current verbosity level

    @return the verbosity level. 
 */
freesasa_verbosity
freesasa_get_verbosity(void);

/**
    Allocate empty structure.

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @return Empty structure. NULL if memory allocation failure.
 */
freesasa_structure*
freesasa_structure_new(void);

/**
    Init structure with coordinates from pdb-file.

    Reads in a PDB-file and generates a structure object.
    Automatically skips hydrogens. If an atom has alternative
    coordinates, only the first alternative is used. If a file has
    more than one `MODEL` (as in NMR structures) only the first model
    is used. User specifies if `HETATM` entries and/or hydrogen atoms
    should be included. It is also possible to specify that all MODELs
    should be joined to one large structure. If more fine-grained
    control over which atoms to include is needed, the PDB-file needs
    to be modified before calling this function, or atoms can be added
    manually one by one using freesasa_structure_add_atom().

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @param pdb Input PDB-file.
    @param options Bitfield. 0 means only use non-hydrogen `ATOM`
      entries, first MODEL only. ::FREESASA_INCLUDE_HETATM and
      ::FREESASA_INCLUDE_HYDROGEN can be used to include more
      atoms. ::FREESASA_JOIN_MODELS can be used if input has several
      models that should be considered part of the same structure. The
      options can be included using `|`, for example
      `FREESASA_INCLUDE_HETATM | FREESASA_INCLUDE_HYDROGEN` means
      include both hydrogens and HETATMs.
    @return The generated structure. Returns `NULL` and prints error
      if input is invalid or  memory allocation failure.
 */

freesasa_structure*
freesasa_structure_from_pdb(FILE *pdb,
                            int options);
/**
    Init array of structures from PDB.
    
    Either iniatilize one structure per model in multimodel PDB, or
    one per chain, or both. Otherwise equivalent to
    freesasa_structure_from_pdb(). 

    Returns dynamically allocated array of size n. Its members should
    be freed using freesasa_structure_free() and the array itself with
    free().

    @param pdb Input PDB-file.  
    
    @param n Number of structures found are written to this integer.
    
    @param options Bitfield. 0 means only use non-hydrogen `ATOM`
      entries. ::FREESASA_INCLUDE_HETATM and
      ::FREESASA_INCLUDE_HYDROGEN can be used to include more
      atoms. ::FREESASA_SEPARATE_MODELS and
      ::FREESASA_SEPARATE_CHAINS can be used to generate one structure
      per model and one structure per chain, respectively.  All four
      options can be combined using `|`, analogously to
      freesasa_structure_from_pdb().

    @return Array of structures. Prints error message(s) and returns
      NULL if there were problems reading input or a memory allocation
      failure.
 */
freesasa_structure**
freesasa_structure_array(FILE *pdb,
                         int *n,
                         int options);

/**
    Add individual atom to structure.

    A structure can be built by adding atoms one by one. Storing
    residue numbers as strings allows for non-numeric labels. Will
    include hydrogens if added (i.e. up to caller to make sure these
    are excluded if necessesary).

    The three string arguments all have specific lengths, specified by
    the corresponding PDB ATOM fields. It is important to use the same
    padding as in the PDB specification (i.e. `" CA "` instead of
    `"CA"` or `" CA"`). (This might be more flexible in future
    versions of the library).

    @param structure The structure to add to.
    @param atom_name String of 4 characters, of the format `" CA "`, `" OXT"`, etc.
    @param residue_name String of 3 charachters, of the format `"ALA"`, `"PHE"`, etc.
    @param residue_number String of 4 characters, of the format `"   1"`, `" 123"`, etc.
    @param chain_label Any character to label chain, typically `'A'`, `'B'`, etc.
    @param x x-coordinate of atom.
    @param y y-coordinate of atom.
    @param z z-coordinate of atom.

    @return ::FREESASA_SUCCESS if input valid. ::FREESASA_FAIL if any
    of the strings are malformatted or if memoray allocation
    fails. ::FREESASA_WARN if the atom type is unknown.
 */
int
freesasa_structure_add_atom(freesasa_structure *structure,
                            const char* atom_name,
                            const char* residue_name,
                            const char* residue_number,
                            char chain_label,
                            double x, double y, double z);

/**
    Create new structure consisting of a selection chains from the provided structure.

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @param structure Input structure.
    @param chains String of chain labels (e.g. "AB")

    @return A new structure consisting only of the specified
    chains. Returns NULL if the input structure doesn't have any
    matching chains or if memory allocation fails.
 */
freesasa_structure*
freesasa_structure_get_chains(const freesasa_structure *structure, 
                              const char* chains);

/**
    Get string listing all chains in structure.

    @param structure The structure.
    @return String with all chain labels in structure ("A", "ABC", etc).
 */
const char*
freesasa_structure_chain_labels(const freesasa_structure *structure);

/**
    Get number of atoms.

    @param structure The structure.
    @return Number of atoms.
 */
int
freesasa_structure_n(const freesasa_structure *structure);

/**
    Free structure.

    @param structure The structure to free.
 */
void
freesasa_structure_free(freesasa_structure* structure);

/**
    Calculates radii of all atoms in the structure using provided
    classifier.

    Return value is dynamically allocated, should be freed with
    standard free().

    @param structure The structure.  
    @param classifier The classifier. If NULL the default is used.  
    @return Array of radii. NULL if there are unrecognized atoms or if
    memory allocation fails.
 */
double*
freesasa_structure_radius(const freesasa_structure *structure,
                          const freesasa_classifier *classifier);

/**
    Get atom name.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Atom name in the form `" CA "`, `" OXT"`, etc.
 */
const char*
freesasa_structure_atom_name(const freesasa_structure *structure,
                             int i);

/**
    Get residue name.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Residue name in the form `"ALA"`, `"PHE"`, etc.
 */
const char*
freesasa_structure_atom_res_name(const freesasa_structure *structure,
                                 int i);

/**
    Get residue number.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Residue name in the form `"   1"`, `" 123"`, etc.
 */
const char*
freesasa_structure_atom_res_number(const freesasa_structure *structure,
                                   int i);

/**
    Get chain label.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Chain label (`'A'`, `'B'`, etc.)
 */
char
freesasa_structure_atom_chain(const freesasa_structure *structure,
                              int i);

/**
    Get model number for structure.

    Useful if structure was generated with freesasa_structure_array().

    @param structure The structure.  
    @return The model number. 0 means no model number has been read.
 */
int
freesasa_structure_model(const freesasa_structure *structure);

#ifdef __cplusplus
}
#endif

#endif /* FREESASA_H */
