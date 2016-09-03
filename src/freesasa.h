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

//! Atoms can be of the classes Apolar, Polar or Unknown.
typedef enum {
    FREESASA_ATOM_APOLAR=0,
    FREESASA_ATOM_POLAR=1,
    FREESASA_ATOM_UNKNOWN=2,
} freesasa_atom_class;

/**
    Options for reading structure from PDB.  To be combined in options
    bitfield in freesasa_structure_from_pdb(),
    freesasa_structure_array() and freesasa_structure_add_atom_wopt(). 
    See documentation for each function for which options are applicable.
 */
enum freesasa_structure_options {
    FREESASA_INCLUDE_HETATM=1, //!< Include HETATM entries
    FREESASA_INCLUDE_HYDROGEN=1<<2, //!< Include hydrogen atoms
    FREESASA_SEPARATE_MODELS=1<<3, //!< Read MODELs as separate structures
    FREESASA_SEPARATE_CHAINS=1<<4, //!< Read separate chains as separate structures
    FREESASA_JOIN_MODELS=1<<5, //!< Read MODELs as part of one big structure
    FREESASA_HALT_AT_UNKNOWN=1<<6, //!< Halt reading when unknown atom is encountered.
    FREESASA_SKIP_UNKNOWN=1<<7, //!< Skip atom when unknown atom is encountered.
    FREESASA_RADIUS_FROM_OCCUPANCY=1<<8, //!< Read atom radius from occupancy field.
};

//! Controls output format, can be combined in options bitfield in freesasa_export_tree()
enum freesasa_output_options {
    FREESASA_OUTPUT_ATOM=1, //!< Output data for atoms, residues, chains and structure
    FREESASA_OUTPUT_RESIDUE=1<<2, //!< Output data for residues, chains and structure
    FREESASA_OUTPUT_CHAIN=1<<3, //!< Output data for chains and structure
    FREESASA_OUTPUT_STRUCTURE=1<<4, //!< Output data only for the whole structure
    FREESASA_RSA=1<<5, //!< Write RSA output (not affected by atom, residue, etc above)
    FREESASA_JSON=1<<6, //!< Write JSON output
    FREESASA_XML=1<<7, //!< Wite XML output
    //! Don't output relative areas, for example if structure has
    //! manually set radii, invalidating reference values
    FREESASA_OUTPUT_SKIP_REL=1<<8,
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

/**
    Struct for structure object.

    The struct includes the coordinates and radius of each atom, and
    its name, residue-name, etc. If it was initiated from a PDB file
    enough info will be stored so that a PDB-file can be printed using
    the original one as template (see freesasa_write_pdb()).
 */
typedef struct freesasa_structure freesasa_structure;

//! Struct to store results of SASA calculation @ingroup API
typedef struct {
    double total; //!< Total SASA in Ångström^2
    double *sasa; //!< SASA of each atom in Ångström^2
    int n_atoms;  //!< Number of atoms
} freesasa_result;

/**
    Struct to store integrated SASA values for either a full structure
    or a subset thereof.

    Use freesasa_classifier_classify_result() to turn a
    ::freesasa_result into a ::freesasa_nodearea. Each
    ::freesasa_structure_node is associated with a
    ::freesasa_nodearea.
 */
typedef struct {
    const char *name;  //!< Name of substructure
    double total;      //!< Total SASA
    double main_chain; //!< Main-chain/Backbone SASA
    double side_chain; //!< Side-chain SASA
    double polar;      //!< Polar SASA
    double apolar;     //!< Apolar SASA
    double unknown;    //!< SASA of unknown class (neither polar nor apolar)
} freesasa_nodearea;

//! Node types
typedef enum {
    FREESASA_NODE_ATOM, //!< Atom node
    FREESASA_NODE_RESIDUE, //!< Residue node
    FREESASA_NODE_CHAIN, //!< Chain node
    FREESASA_NODE_STRUCTURE, //!< Structure node
    FREESASA_NODE_ROOT, //!< Root node, wraps one or more structures
    FREESASA_NODE_NONE //!< for specifying not a valid node
} freesasa_node_type;

/**
    A node representing a structure, chain, residue or atom in a
    structure (see ::freesasa_node_type). These can be linked in a
    tree created by freesasa_structure_tree(), freed by
    freesasa_structure_tree_free().
 */
typedef struct freesasa_structure_node freesasa_structure_node;


/**
    Struct that can be used to determine classes (polar/apolar) and
    radii of atoms. Initiated from
    freesasa_classifier_from_filename(). The classifiers
    ::freesasa_default_classifier, ::freesasa_protor_classifier,
    ::freesasa_naccess_classifier and ::freesasa_classifier are const
    classifiers that can be used directly.
 */
typedef struct freesasa_classifier freesasa_classifier;

//! Classifier using ProtOr radii and classes
extern const freesasa_classifier freesasa_protor_classifier;

//! Classifier using NACCESS radii and classes
extern const freesasa_classifier freesasa_naccess_classifier;

//! Classifier using OONS radii and classes
extern const freesasa_classifier freesasa_oons_classifier;

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
    Use a classifier to determine the radius of a given atom.

    @param classifier The classifier.
    @param res_name The residue name (ALA/VAL/U/C/...)
    @param atom_name The atom name (CA/N/CB/...)
    @return The radius, negative if atom unknown.
 */
double
freesasa_classifier_radius(const freesasa_classifier *classifier,
                           const char *res_name,
                           const char *atom_name);

/**
    Use a classifier to determine the class of a given atom.

    @param classifier The classifier.
    @param res_name The residue name (ALA/VAL/U/C/...)
    @param atom_name The atom name (CA/N/CB/...)
    @return The class.
 */
freesasa_atom_class
freesasa_classifier_class(const freesasa_classifier *classifier,
                          const char *res_name, 
                          const char *atom_name);

/**
    Names for ::freesasa_atom_class.

    @param atom_class The class.
    @return Name of class.
 */
const char*
freesasa_classifier_class2str(freesasa_atom_class atom_class);

/**
    The name of a classifier.

    @param classifier The classifier.
    @return The name of the classifier.
 */
const char*
freesasa_classifier_name(const freesasa_classifier *classifier);

/**
    Classify results.

    Adds up the SASA of Polar/Apolar/Unknown atoms, and
    main-chain/side-chain atoms for the whole protein. Use
    freesasa_result2tree() for a more fine-grained analysis.

    @param classifier The classifier to use
    @param structure The structure the results are based on
    @param result The results
    @return A struct with all the results.
 */
freesasa_nodearea
freesasa_classifier_classify_result(const freesasa_classifier *classifier,
                                    const freesasa_structure *structure,
                                    const freesasa_result *result);

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
                          const freesasa_result *result,
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
    Write results of claculation to file.

    @param log Output-file.
    @param result SASA values.
    @param name Name of the protein, if NULL "unknown" used.
    @param chains The chains used in the calculation, can be NULL
    @param class_area The SASA values for each class, if NULL
      only total SASA printed
    @return ::FREESASA_SUCCESS on success, ::FREESASA_FAIL if problems
      writing to file.
 */
int
freesasa_write_result(FILE *log,
                      freesasa_result *result,
                      const char *name,
                      const char *chains,
                      const freesasa_nodearea *class_area);
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

      - `options & ::FREESASA_RADIUS_FROM_OCCUPANCY == 1`: Read atomic
         radii from Occupancy field in PDB file.

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

    @param structure A structure.
    @return Number of residues.
 */
int
freesasa_structure_n_residues(const freesasa_structure *structure);

/**
    Get number of chains.

    @param structure A structure.
    @return The number of chains in the structure.
 */
int
freesasa_structure_n_chains(const freesasa_structure *structure);

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
freesasa_structure_set_radius(freesasa_structure *structure,
                              const double* radii);

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
    Get atom radius.

    Asserts that index i is within bounds. 

    @param structure The structure.
    @param i Atom index.
    @return Atom radius.
 */
double
freesasa_structure_atom_radius(const freesasa_structure *structure,
                               int i);
/**
    Set atom radius.
    
    Asserts that index i is within bounds. 

    @param structure The structure.
    @param radius The radius.
    @param i Atom index.
 */
void
freesasa_structure_atom_set_radius(freesasa_structure *structure,
                                   int i,
                                   double radius);

/**
    Get name of residue.

    @param structure The structure.
    @param r_i Residue index (in whole structure)
    @return Name of residue
 */
const char*
freesasa_structure_residue_name(const freesasa_structure *structure,
                                int r_i);

/**
    Get residue number.

    @param structure The structure.
    @param r_i Residue index (in whole structure).
    @return Residue number as string.
 */
const char*
freesasa_structure_residue_number(const freesasa_structure *structure,
                                  int r_i);
/**
    Get chain residue belongs to.

    @param structure The structure.
    @param r_i Residue index (in whole structure).
    @return Chain label.
 */
char
freesasa_structure_residue_chain(const freesasa_structure *structure,
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
 
o    @param structure A structure.
    @param r_i Residue index.
    @param first First atom of residue `r_i` will be stored here.
    @param last Last atom of residue `r_i` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if index `r_i` is invalid.
 */
int
freesasa_structure_residue_atoms(const freesasa_structure *structure,
                                 int r_i, 
                                 int *first,
                                 int *last);

/**
    Get indices of first and last atoms of a chain
 
    @param structure A structure.
    @param chain The chain label.
    @param first First atom of `chain` will be stored here.
    @param last Last atom of `chain` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if `chain` not found.
 */
int
freesasa_structure_chain_atoms(const freesasa_structure *structure,
                               char chain,
                               int *first,
                               int *last);

/**
    Get indices of first and last residues of a chain

    @param structure A structure.
    @param chain The chain label.
    @param first First residue of `chain` will be stored here.
    @param last Last residue of `chain` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if `chain` not found.
 */
int
freesasa_structure_chain_residues(const freesasa_structure *structure,
                                  char chain,
                                  int *first,
                                  int *last);

/**
    Generates a tree that represents the structure, with the levels
    described by ::freesasa_node_type. 

    Use freesasa_structure_node_children() to traverse the tree, and
    freesasa_structure_node_type(), freesasa_structure_node_area(),
    and freesasa_node_name() to explore the properties of each node.

    The tree should be freed using freesasa_structure_tree_free().

    @param result SASA values for the structure
    @param structure The structure the results are based.
    @param classifier Classifier to determine which atoms are polar
      and apolar, and, if available, give reference values to
      calculate relative SASAs for residues. If NULL
      ::freesasa_default_classifier will be used. For consistent 
      results this should be the same classifier used to determine
      atomic radii when the structure was created.
    @param name The name of the structure
    @return The root node of the tree. NULL if memory allocation fails.
 */
freesasa_structure_node *
freesasa_result2tree(const freesasa_result *result,
                     const freesasa_structure *structure,
                     const freesasa_classifier *classifier,
                     const char *name);

/**
    Join two ::structure_node-trees.

    Allows joining several calculations into one output file.

    @param tree1 The joint tree will be stored here
    @param tree2 Will be added to tree1, and then changed to NULL,
      since ownership of its contents have been transferred to tree1.
    @return ::FREESASA_SUCCESS upon sucecss. ::FREESASA_WARN if tree1
      and tree2 were created using different classifiers.
 */
int
freesasa_structure_node_join_trees(freesasa_structure_node *tree1,
                                   freesasa_structure_node **tree2);

/**
    Outputs result in format specified by options.

    @param output Output file.
    @param root Structure tree containing results (generated using 
      freesasa_result2tree()).
    @param parameters Parameters used in the calculated, printed for 
      reference in some formats.
    @param options Bitfield specifying output format, see 
      ::freesasa_output_options.
*/
int
freesasa_export_tree(FILE *output,
                     const freesasa_structure_node *root,
                     const freesasa_parameters *parameters,
                     int options);

/**
    Free ::freesasa_structure_node-tree.

    Will not free anything if the node has a parent, i.e. if this node
    is the not the root of the tree it belongs too.

    @param root Root of the tree to free
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if the node hasa a
      parent.
 */
int
freesasa_structure_node_free(freesasa_structure_node *root);

/**
    The ::freesasa_nodearea of all atoms belonging to a node.

    Only works if the areas have been added using
    freesasa_structure_tree_fill().

    @param node The node. 
    @return The area. NULL if no area has been attached to this node.
 */
const freesasa_nodearea *
freesasa_structure_node_area(const freesasa_structure_node *node);

/**
    The children of a node.

    Use freesasa_structure_node_next() to access next sibling.

    @param node The node.  
    @return Pointer to the first child of a node. NULL if the node has no
      children.
 */
const freesasa_structure_node *
freesasa_structure_node_children(const freesasa_structure_node *node);

/**
    Next sibling of a node.

    @param node The node.
    @return The next node, NULL if this is the last node.
 */
const freesasa_structure_node *
freesasa_structure_node_next(const freesasa_structure_node *node);

/**
    The parent of a node.

    @param node The node.
    @return The parent node. NULL if the node has no parent.
 */
const freesasa_structure_node *
freesasa_structure_node_parent(const freesasa_structure_node *node);

/**
    The type of a node.

    @param node The node.
    @return The ::freesasa_node_type.
 */
freesasa_node_type
freesasa_structure_node_type(const freesasa_structure_node *node);

/**
    The name of a node.
    
    @param node The node.
    @return The name. NULL if the node has no name.
 */
const char *
freesasa_structure_node_name(const freesasa_structure_node *node);

/**
    The name of the classifier used to generate the node.

    @param node The node (has to be of type ::FREESASA_NODE_ROOT)
    @return The name of the classifier
 */
const char*
freesasa_structure_node_classified_by(const freesasa_structure_node *node);

/**
    Is atom polar.

    @param node The atom (has to be of type ::FREESASA_NODE_ATOM).
    @return 1 if polar, 0 else.
 */
int
freesasa_structure_node_atom_is_polar(const freesasa_structure_node *node);

/**
    Does atom belong to the main chain/backbone.

    @param node The atom (has to be of type ::FREESASA_NODE_ATOM).
    @return 1 if mainchain, 0 else.
 */
int
freesasa_structure_node_atom_is_mainchain(const freesasa_structure_node *node);

/**
    Atom radius.

    @param node The atom (has to be of type ::FREESASA_NODE_ATOM).
    @return The radius.
 */
double
freesasa_structure_node_atom_radius(const freesasa_structure_node *node);

/**
    Residue number.

    @param node The residue (has to be of type ::FREESASA_NODE_RESIDUE).
    @return String with residue number.
 */
const char *
freesasa_structure_node_residue_number(const freesasa_structure_node *node);

/**
    Number of atoms in a residue.

    @param node The residue (has to be of type ::FREESASA_NODE_RESIDUE).
    @return Number of atoms.
 */
int
freesasa_structure_node_residue_n_atoms(const freesasa_structure_node *node);

/**
    The reference area for a node from the classifier used to
    generate the tree.

    @param node The node (has to be of type ::FREESASA_NODE_RESIDUE)
    @return The reference area. NULL if area not available or if node
      is not a residue.
 */
const freesasa_nodearea *
freesasa_structure_node_residue_reference(const freesasa_structure_node *node);

/**
    The number of residues in a chain.

    @param node The chain (has to be of type ::FREESASA_NODE_CHAIN).
    @return Number of residues.
 */
int
freesasa_structure_node_chain_n_residues(const freesasa_structure_node *node);

/**
    The number of chains in a structure.

    @param node The structure (has to be of type ::FREESASA_NODE_STRUCTURE).
    @return Number of chains.
 */
int
freesasa_structure_node_structure_n_chains(const freesasa_structure_node *node);

/**
    All chain labels in a structure.

    @param node The structure (has to be of type ::FREESASA_NODE_STRUCTURE).
    @return Chain labels as null-terminated string.
 */
const char *
freesasa_structure_node_structure_chain_labels(const freesasa_structure_node *node);

// Deprecated functions below, from 1.x API

/**
    Struct used to store n string-value-pairs (strvp) in arrays of
    doubles and strings. freesasa_strvp_free() assumes both arrays
    and strings are dynamically allocated.

    @deprecated use freesasa_nodearea instead
 */
typedef struct {
    double *value; //!< Array of values
    char **string; //!< Array of strings
    int n;         //!< Number of values and strings
} freesasa_strvp;

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

    @deprecated Use freesasa_classifier_classify_result() instead.
 */
freesasa_strvp*
freesasa_result_classify(const freesasa_result *result,
                         const freesasa_structure *structure,
                         const freesasa_classifier *classifier);


/**
    Frees a ::freesasa_strvp object

    @param strvp the object to be freed

    @deprecated
 */
void
freesasa_strvp_free(freesasa_strvp *strvp);

#ifdef __cplusplus
}
#endif

#endif /* FREESASA_H */
