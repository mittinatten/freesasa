#ifndef FREESASA_H
#define FREESASA_H

/**
    @file
    @author Simon Mitternacht
    @copyright [MIT License](md_license.html)

    @brief Functions and datatypes for performing and analyzing SASA
    calculations.

    This header provides the functions and data types necessary to
    perform and analyze a SASA calculation using FreeSASA, including
    facilities to customize assignment of radii to, and classification
    of, atoms. There are also functions to access properties of a
    structure, to allow refined analysis of the results. The page @ref
    API shows how to set up and perform a simple SASA calculation.

 */

#include <stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

//! The FreeSASA algorithms. 
typedef enum {
    FREESASA_LEE_RICHARDS, //!< Lee & Richards' algorithm
    FREESASA_SHRAKE_RUPLEY //!< Shrake & Rupley's algorithm
} freesasa_algorithm;

//! Verbosity levels. @see freesasa_set_verbosity() @see freesasa_get_verbosity()
typedef enum {
    FREESASA_V_NORMAL, //!< Print all errors and warnings.
    FREESASA_V_NOWARNINGS, //!< Print only errors.
    FREESASA_V_SILENT, //!< Print no errors and warnings.
    FREESASA_V_DEBUG, //!< Print all errors, warnings and debug messages.
} freesasa_verbosity;

// Default parameters
#define FREESASA_DEF_ALGORITHM FREESASA_LEE_RICHARDS //!< Default algorithm
#define FREESASA_DEF_PROBE_RADIUS 1.4 //!< Default probe radius (in Ångström).
#define FREESASA_DEF_SR_N 100 //!< Default number of test points in S&R.
#define FREESASA_DEF_LR_N 20 //!< Default number of slices per atom  in L&R.

//! Default ::freesasa_classifier
#define freesasa_default_classifier freesasa_protor_classifier

//! Default ::freesasa_rsa_reference
#define freesasa_default_rsa freesasa_protor_rsa

//! Default number of threads. Value will depend on if library was
//! compiled with or without thread support. (2 with threads, 1
//! without)
const extern int FREESASA_DEF_NUMBER_THREADS; 

//! Error codes. Can rely upon FREESASA_SUCCESS being 0 and the errors
//! having negative numbers.
enum freesasa_error_codes {
    FREESASA_SUCCESS=0, //!< All is ok (value will always be zero).
    FREESASA_FAIL=-1, //!< Something went seriously wrong (value will always be negative).
    FREESASA_WARN=-2, //!< Something went wrong, but results might still be meaningful (value will always be negative).
};

/**
    Options for reading structure from PDB.  To be combined in options
    bitfield in freesasa_structure_from_pdb(),
    freesasa_structure_array() and freesasa_structure_add_atom_wopt(). 
    See documentation for each function for which options are applicable.
 */
enum freesasa_structure_options {
    FREESASA_INCLUDE_HETATM=1, //!< Include HETATM entries
    FREESASA_INCLUDE_HYDROGEN=2, //!< Include hydrogen atoms
    FREESASA_SEPARATE_MODELS=4, //!< Read MODELs as separate structures
    FREESASA_SEPARATE_CHAINS=8, //!< Read separate chains as separate structures
    FREESASA_JOIN_MODELS=16, //!< Read MODELs as part of one big structure
    FREESASA_HALT_AT_UNKNOWN=32, //!< Halt reading when unknown atom is encountered.
    FREESASA_SKIP_UNKNOWN=64, //!< Skip atom when unknown atom is encountered.
};

//! The maximum length of a selection name @see freesasa_select_area()
#define FREESASA_MAX_SELECTION_NAME 50

//! Struct to store parameters for SASA calculation @ingroup API
typedef struct {
    freesasa_algorithm alg;       //!< Algorithm
    double probe_radius;          //!< Probe radius (in Ångström)
    int shrake_rupley_n_points;   //!< Number of test points in S&R calculation
    int lee_richards_n_slices;    //!< Number of slices per atom in L&R calculation
    int n_threads;                //!< Number of threads to use, if compiled with thread-support
} freesasa_parameters;

//! The default parameters for FreeSASA
extern const freesasa_parameters freesasa_default_parameters;

//! Struct to store results of SASA calculation @ingroup API
typedef struct {
    double total; //!< Total SASA in Ångström^2
    double *sasa; //!< SASA of each atom in Ångström^2
    int n_atoms;  //!< Number of atoms
} freesasa_result;

//! Struct to store SASA values for a named residue
typedef struct {
    const char *name;  //!< Residue name
    double total;      //!< Total SASA
    double main_chain; //!< Main-chain/Backbone SASA
    double side_chain; //!< Side-chain SASA
    double polar;      //!< Polar SASA
    double apolar;     //!< Apolar SASA
} freesasa_residue_sasa;

/**
    Struct for structure object.

    The struct includes the coordinates and radius of each atom, and
    its name, residue-name, etc. If it was initiated from a PDB file
    enough info will be stored so that a PDB-file can be printed using
    the original one as template (see freesasa_write_pdb()).
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

    /**
        Function that returns an atom radius. Should return negative
        value if atom not recognized.
    */
    double (*radius)(const char* res_name,
                     const char* atom_name,
                     const struct freesasa_classifier*);

    /**
        Function that returns the class [0,1,...,n_classes-1] of an
        atom, should return ::FREESASA_WARN if atom not recognized.
    */
    int (*sasa_class)(const char* res_name,
                      const char* atom_name,
                      const struct freesasa_classifier*);

    //! Function that converts a class to its string descriptor.
    const char* (*class2str)(int the_class,
                             const struct freesasa_classifier*);

    //! Function that can be called to free the config-pointer
    void (*free_config)(void*);
} freesasa_classifier;

//! Classifier using ProtOr radii and classes
extern const freesasa_classifier freesasa_protor_classifier;

//! Classifier using NACCESS radii and classes
extern const freesasa_classifier freesasa_naccess_classifier;

//! Classifier using OONS radii and classes
extern const freesasa_classifier freesasa_oons_classifier;

//! Used as reference in generation of RSA file
typedef struct {
    //! Name of RSA reference
    char *name;

    /**
        Array containing max SASA-values for named residues types.
        Last element of array should have name == NULL.
    */
    const freesasa_residue_sasa *max;

    /**
        A classifier whose sasa_class() function returns 0 for apolar
        atoms and non-zero for polar atoms. This is true for
        ::freesasa_protor_classifier, ::freesasa_oons_classifier and
        ::freesasa_naccess_classifier.
    */
    const freesasa_classifier *polar_classifier;
    
    /**
       A classifier whose sasa_class() function returns 0 for side
       chain atoms and non-zero for backbone.
    */
    const freesasa_classifier *bb_classifier;
} freesasa_rsa_reference;

//! An RSA-reference for the ProtOr configuration
extern const freesasa_rsa_reference freesasa_protor_rsa;

//! An RSA-reference for the NACCESS configuration
extern const freesasa_rsa_reference freesasa_naccess_rsa;

/**
    Calculates SASA based on a given structure.

    Return value is dynamically allocated, should be freed with
    freesasa_result_free().

    @param structure The structure
    @param parameters Parameters for the calculation, if NULL
    defaults are used.

    @return The result of the calculation, NULL if something went wrong.
 */
freesasa_result *
freesasa_calc_structure(const freesasa_structure *structure,
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
freesasa_result *
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
    Frees a classifier object

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
    Get area of a selection.

    Uses subset of the select syntax from Pymol (name, symbol, resn,
    resi and chain), the keyword "select" is implicit. All commands
    are case insensitive. Valid selections would be, for example,

        selection_name, resn ala+arg 
        selection_name, chain a and resi 1+3-20 and not resn gly

    After selecting the atoms from the ::freesasa_structure pointer
    specified by the command the area of those atoms is summed up
    using the ::freesasa_result pointer.

    @see @ref Selection

    @param command The selection
    @param name The name of the selection is stored here, it should be
      able to store a string of length ::FREESASA_MAX_SELECTION_NAME.
    @param area The area of the selection is stored here
    @param structure The structure to select from
    @param result The results to integrate

    @return ::FREESASA_SUCCESS upon successful selection.
       ::FREESASA_WARN if some illegal selections that could be
       ignored were encountered (see printed
       warnings). ::FREESASA_FAIL if syntax error or memory failure.
 */
int
freesasa_select_area(const char *command,
                     char *name,
                     double *area,
                     const freesasa_structure *structure,
                     const freesasa_result *result);

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
    using freesasa_structure_from_pdb() or freesasa_structure_array().

    @param output File to write to.
    @param result SASA values.
    @param structure Structure to use to print PDB.

    @return ::FREESASA_SUCCESS if file written successfully.
    ::FREESASA_FAIL if there is no previous PDB input to base output
    on or if there were problems writing to the file.
 */
int
freesasa_write_pdb(FILE *output, 
                   freesasa_result *result,
                   const freesasa_structure *structure);

/**
    Print SASA for each chain.

    Prints the contribution of each chain to the total SASA. Each line
    in the output is prefixed by the string `CHAIN`.

    @param output Output file.
    @param result SASA values.
    @param structure The structure. 
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.
 */
int
freesasa_per_chain(FILE *output,
                   freesasa_result *result,
                   const freesasa_structure *structure);

/**
    Print SASA for all residue types to file.

    Prints name/value-pairs with the total SASA of each residue
    type. The standard 20 amino acids are always included in output,
    non-standard ones and nucleotides only if they were present in
    input. Each line in the output is prefixed by the string "`RES`".

    @param output Output file.
    @param result SASA values.
    @param structure The structure.
    @return ::FREESASA_FAIL if problems writing to
      `output`. ::FREESASA_SUCCESS else.
 */
int
freesasa_per_residue_type(FILE *output,
                          freesasa_result *result,
                          const freesasa_structure *structure);

/**
    Print SASA for each residue in the sequence to file.

    Each line in the output is prefixed by the string "`SEQ`".

    @param output Output file.
    @param result SASA values.
    @param structure The structure.
    @return ::FREESASA_FAIL if problems writing to
      `output`. ::FREESASA_SUCCESS else.
 */
int
freesasa_per_residue(FILE *output,
                     freesasa_result *result,
                     const freesasa_structure *structure);

/**
    Log calculation results.

    Prints log of calculation to specified file, equivalent of calling
    first freesasa_write_parameters() and then freesasa_write_result()

    @param log Output-file.
    @param result SASA values.
    @param parameters Parameters to print, if NULL defaults are used
    @param name Name of the protein, if NULL "unknown" used.
    @param class_sasa The SASA values for each class, if NULL
      only total SASA printed
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if problems
      writing to file.

    @deprecated Use freesasa_write_parameters() and
    freesasa_write_result() instead.
 */
int
freesasa_log(FILE *log,
             freesasa_result *result,
             const char *name,
             const freesasa_parameters *parameters,
             const freesasa_strvp* class_sasa);
/**
    Write results of claculation to file.

    @param log Output-file.
    @param result SASA values.
    @param name Name of the protein, if NULL "unknown" used.
    @param chains The chains used in the calculation, can be NULL
    @param class_sasa The SASA values for each class, if NULL
      only total SASA printed
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if problems
      writing to file.
 */
int
freesasa_write_result(FILE *log,
                      freesasa_result *result,
                      const char *name,
                      const char *chains,
                      const freesasa_strvp* class_sasa);
/**
    Print parameters to file

    @param log Output-file
    @param parameters Parameters to print, if NULL defaults are used
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if problems
      writing to file.
*/
int freesasa_write_parameters(FILE *log,
                              const freesasa_parameters *parameters);

/**
    Print RSA-file

    Uses reference SASA values calculated using the default
    configuration (ProtOr radii; Carbon is apolar, all other elements
    polar; probe radius = 1.4 Å). The Ala-X-Ala configurations
    supplied in the directory `rsa` were used as input (the reference
    values themselves are stored statically in the code). At the
    moment there is no support for outputting RSA files for other
    configurations.

    @remark This is still an experimental feature, and the interface
      may still be subject to change without warning.

    @param output Output-file
    @param result SASA values
    @param structure The structure
    @param name Name of the protein
    @param reference Reference to calculate RSA from, if NULL defaults
      are used.
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if problems
      writing to file, or if `structure` is inconsistent.
 */
int
freesasa_write_rsa(FILE *output,
                   const freesasa_result *result,
                   const freesasa_structure *structure,
                   const char *name,
                   const freesasa_rsa_reference *reference);
/**
    Set the global verbosity level.

    @param v the verbosity level
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
    Set where to write errors.
    
    By default stderr is used, this function can be called to redirect
    error output elsewhere.

    @param err The file to write to. If NULL stderr will be used.
 */
void
freesasa_set_err_out(FILE *err);

/**
    Get pointer to error file.

    NULL means stderr is used.
  
    @return The error file.
 */
FILE *
freesasa_get_err_out();

/**
    Allocate empty structure.

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @return Empty structure. NULL if memory allocation failure.
 */
freesasa_structure*
freesasa_structure_new(void);

/**
    Free structure.

    @param structure The structure to free.
 */
void
freesasa_structure_free(freesasa_structure* structure);

/**
    Init structure with coordinates from pdb-file.

    Reads in a PDB-file and generates a structure object.
    Automatically skips hydrogens and HETATM. If an atom has
    alternative coordinates, only the first alternative is used. If a
    file has more than one `MODEL` (as in NMR structures) only the
    first model is used. The provided classifier is used to determine
    the radius of each atom, if the atom is not recognized the element
    of the atom is guessed, and that element's VdW radius used. If
    this fails its radius is set 0, which means that it won't
    contribute to SASA, but a radius from another source can be set
    through freesasa_structure_set_radius(). All these behaviors can
    be modified through the `options` bitfield argument:

      - `options == 0`:
         Default behavior

      - `options & ::FREESASA_INCLUDE_HYDROGEN == 1`:
         Include hydrogens.

      - `options & ::FREESASA_INCLUDE_HETATM == 1`:
         Include HETATM.

      - `options & ::FREESASA_JOIN_MODELS == 1`:
         Join models.

      - `options & ::FREESASA_SKIP_UNKNOWN == 1`: Skip unknown
         atoms.
      
      - `options & ::FREESASA_HALT_AT_UNKNOWN == 1`: Halt at unknown
         atom and return NULL. Overrides ::FREESASA_SKIP_UNKNOWN.

    If a more fine-grained control over which atoms to include is
    needed, the PDB-file needs to be modified before calling this
    function, or atoms can be added manually one by one using
    freesasa_structure_add_atom() or
    freesasa_structure_add_atom_wopt().

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @param pdb A PDB file

    @param classifier A freesasa_classifier to determine radius of
      atom. If NULL default classifier is used.

    @param options A bitfield to determine what atoms to include and what to do
      when atoms are not recognized by classifier.
      
    @return The generated structure. Returns `NULL` and prints error
      if input is invalid or  memory allocation failure.
 */

freesasa_structure*
freesasa_structure_from_pdb(FILE *pdb,
                            const freesasa_classifier *classifier,
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
    @param classifier A classifier to calculate atomic radii.
    @param options Bitfield. Either or both of
      ::FREESASA_SEPARATE_MODELS and ::FREESASA_SEPARATE_CHAINS can be
      used to generate one structure per model and one structure per
      chain, respectively (will return NULL if neither is
      specified). See freesasa_structure_from_pdb() for documentation
      on options for deciding what atoms to include
      (::FREESASA_JOIN_MODELS is not supported here).
    @return Array of structures. Prints error message(s) and returns
      NULL if there were problems reading input, if invalid value of
      `options`, or upon a memory allocation
      failure.
 */
freesasa_structure **
freesasa_structure_array(FILE *pdb,
                         int *n,
                         const freesasa_classifier *classifier,
                         int options);

/**
    Add individual atom to structure using default behavior.
    
    Equivalent to calling freesasa_structure_add_atom_wopt(), with
    `classifier = NULL` and `options = 0`.
    
    @param structure The structure to add to.
    @param atom_name String of 4 characters, of the format `" CA "`, `" OXT"`, etc.
    @param residue_name String of 3 charachters, of the format `"ALA"`, `"PHE"`, etc.
    @param residue_number String of 4 characters, of the format `"   1"`, `" 123"`, etc.
    @param chain_label Any character to label chain, typically `'A'`, `'B'`, etc.
    @param x x-coordinate of atom.
    @param y y-coordinate of atom.
    @param z z-coordinate of atom.
    @return ::FREESASA_SUCCESS on normal execution. ::FREESASA_FAIL if
      if memory allocation fails.
 */
int 
freesasa_structure_add_atom(freesasa_structure *structure,
                            const char* atom_name,
                            const char* residue_name,
                            const char* residue_number,
                            char chain_label,
                            double x, double y, double z);
/**
    Add individual atom to structure.

    A structure can be built by adding atoms one by one. Storing
    residue numbers as strings allows for non-numeric labels. Will
    include hydrogens if added (i.e. up to caller to make sure these
    are excluded if necessesary).

    The atom name, residue name, etc are checked by the classifier,
    and depending on the value of `options` different things will
    happen when unknown atoms are encountered. In all cases the user
    will be alerted of what has happened through warnings or error
    messages:

      - `options == 0` means guess element of unknown atoms, and use
        that element's VdW radius. If this fails, assign radius 0. A 0
        radius means this atom won't contribute to the SASA, but will
        still be there if we want to use
        freesasa_structure_set_radius() to assign a radius from
        another source.

      - `options & FREESASA_SKIP_UNKNOWN == 1` skip unknown atoms,
        return ::FREESASA_WARN.

      - `options & FREESASA_HALT_AT_UNKNOWN == 1` skip unknown atoms,
        return ::FREESASA_FAIL. Overrides ::FREESASA_SKIP_UNKNOWN.

    @see Because the argument list is so long, freesasa_structure_add_atom()
         is a shortcut to call this with defaults.

    @param structure The structure to add to.
    @param atom_name The atom name: `" CA "`,`"CA"`, `" OXT"`, etc.
    @param residue_name The residue name: `"ALA"`, `"PHE"`, etc.
    @param residue_number String of 4 characters, of the format `"   1"`, `" 123"`, etc.
    @param chain_label Any character to label chain, typically `'A'`, `'B'`, etc.
    @param x x-coordinate of atom.
    @param y y-coordinate of atom.
    @param z z-coordinate of atom.
    @param classifier A freesasa_classifier to determine radius of atom and to
      decide if to keep atom or not (see options).
    @param options A bitfield to determine what to do with unknown atoms (see above).
      
    @return ::FREESASA_SUCCESS on normal execution. ::FREESASA_FAIL if
       if memory allocation fails or if halting at unknown
       atom. ::FREESASA_WARN if skipping atom.
 */
int 
freesasa_structure_add_atom_wopt(freesasa_structure *structure,
                                 const char *atom_name,
                                 const char *residue_name,
                                 const char *residue_number,
                                 char chain_label,
                                 double x, double y, double z,
                                 const freesasa_classifier *classifier,
                                 int options);

/**
    Create new structure consisting of a selection chains from the
    provided structure.

    Simply looks for chain labels that match the characters in the
    provided string.

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
    Get number of residues.

    Calculated crudely by determining the number of unique
    combinations of residue number and chain label contained in the
    structure. If residues are mingled i.e. atoms of the same residue
    are in non-contiguous regions of the file, this function might be
    off.

    @param s A structure.
    @return Number of residues.
 */
int
freesasa_structure_n_residues(const freesasa_structure *s);

/**
    Get number of chains.

    @param s A structure.
    @return The number of chains in the structure.
 */
int
freesasa_structure_n_chains(const freesasa_structure *s);

/**
    Returns a pointer to an array of the radii of each atom.

    @param structure The structure.  
    @return Array of radii. If NULL structure has not been properly
      initialized.
 */
const double*
freesasa_structure_radius(const freesasa_structure *structure);

/**
    Override the radii set when the structure was initialized.

    Makes a copy of the provided array.

    @param structure The structure.
    @param radii An array of radii, should have same dimension
      as the number of atoms in the structure.
 */
void
freesasa_structure_set_radius(freesasa_structure *structure, const double* radii);

/**
    Get atom name.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.

    @return Atom name in the form `" CA "`, `" OXT"`, etc, if
       structure was initialized from a PDB file, or in whatever form
       it was added through freesasa_structure_add_atom() or
       freesasa_structure_add_atom_wopt().
 */
const char*
freesasa_structure_atom_name(const freesasa_structure *structure,
                             int i);

/**
    Get residue name.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Residue name in the form `"ALA"`, `"PHE"`, etc, if
       structure was initialized from a PDB file, or in whatever form
       it was added through freesasa_structure_add_atom() or
       freesasa_structure_add_atom_wopt().
 */
const char*
freesasa_structure_atom_res_name(const freesasa_structure *structure,
                                 int i);

/**
    Get residue number.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Residue name in the form `"   1"`, `" 123"`, etc, if
       structure was initialized from a PDB file, or in whatever form
       it was added through freesasa_structure_add_atom() or
       freesasa_structure_add_atom_wopt().
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
    Get atom symbol.

    If the structure was initialized from a PDB file the symbol field
    of that file is used. Otherwise the symbol is guessed from atom and
    residue name.

    Asserts that index i is within bounds. 

    @param structure The structure.
    @param i Atom index.
    @return Atom symbol (`" C"`, `" N"`, `"SE"`,etc); 
 */
const char*
freesasa_structure_atom_symbol(const freesasa_structure *structure,
                               int i);

/**
    Get name of residue.

    @param s The structure.
    @param r_i Residue index (in whole structure)
    @return Name of residue
 */
const char*
freesasa_structure_residue_name(const freesasa_structure *s,
                                int r_i);

/**
    Get residue number.

    @param s The structure.
    @param r_i Residue index (in whole structure).
    @return Residue number as string.
 */
const char*
freesasa_structure_residue_number(const freesasa_structure *s,
                                  int r_i);
/**
    Get chain residue belongs to.

    @param s The structure.
    @param r_i Residue index (in whole structure).
    @return Chain label.
 */
char
freesasa_structure_residue_chain(const freesasa_structure *s,
                                 int r_i);

/**
    Get model number for structure.

    Useful if structure was generated with freesasa_structure_array().

    @param structure The structure.  
    @return The model number. 0 means no model number has been read.
 */
int
freesasa_structure_model(const freesasa_structure *structure);

/**
    Get array of coordinates.

    Size of array is 3*N, order of coordinates `x1, y1, z1, ...`.

    @param structure The structure.
    @return Array of coordinates. NULL if structure empty. Size can be
      accessed through freesasa_structure_n() (multiply by three).
 */
const double*
freesasa_structure_coord_array(const freesasa_structure *structure);

/**
    Get indices of first and last atoms of a residue
 
    @param s A structure.
    @param r_i Residue index.
    @param first First atom of residue `r_i` will be stored here.
    @param last Last atom of residue `r_i` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if index `r_i` is invalid.
 */
int
freesasa_structure_residue_atoms(const freesasa_structure *s,
                                 int r_i, 
                                 int *first,
                                 int *last);

/**
    Get indices of first and last atoms of a chain
 
    @param s A structure.
    @param chain The chain label.
    @param first First atom of `chain` will be stored here.
    @param last Last atom of `chain` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if `chain` not found.
 */
int
freesasa_structure_chain_atoms(const freesasa_structure *s,
                                 char chain,
                                 int *first,
                                 int *last);


#ifdef __cplusplus
}
#endif

#endif /* FREESASA_H */
