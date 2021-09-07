#ifndef FREESASA_H
#define FREESASA_H

/**
    @file
    @author Simon Mitternacht
    @copyright [MIT License](md_license.html)

    @brief Functions and datatypes for performing and analyzing SASA
    calculations.

    This header provides the functions and data types necessary to
    perform and analyze a SASA calculation using FreeSASA. They are
    all listed on this page, but also divided into @ref core, @ref
    structure, @ref classifier, @ref node and @ref selection
    modules. The page @ref API shows how to set up and perform a
    simple SASA calculation.

    @defgroup core Core

    The core functions to perform a calculation

    @defgroup node Node

    Represent results as a tree.

    Results are represented hierarchically as a tree with the
    following levels (see also ::freesasa_node_type)

    - Root: wrapper, to allow joining results of several calculations.
    - Result: wrapper for an individual result, contains information
      about calculation parameters, the classifier used and the input.
    - Structure: an individual molecule, the highest level where
      actual SASA values are stored.
    - Chain: an individual chain.
    - Residue: typically an amino acid or nucleic acid.
    - Atom: lowest level.

    The tree can be traversed using freesasa_node_children(),
    freesasa_node_next() and freesasa_node_parent(). The type of a
    node is determined by freesasa_node_type(). There are some
    properties that are common to all or most levels of the
    tree. These accessors simply have the prefix `freesasa_node`, like
    for example freesasa_node_name() and the above-mentioned
    functions. Other properties are specific to a special level, and
    then the prefix of the accessor functions will be
    `freesasa_node_atom` or `freesasa_node_structure`, etc. The class
    of an atom can for example be accessed using
    freesasa_node_atom_is_polar().

    The nodes of a tree are to be considered read-only, and changes
    are made only to the root node, initialized using
    freesasa_tree_new() or freesasa_tree_init(), and modified using
    freesasa_tree_add_result(), freesasa_tree_join(). The one
    exception where a lower level node can be modified is
    freesasa_node_structure_add_selection().

    @defgroup structure Structure

    @brief Representation of macromolecular structures.

    Interface for macromolecule structures, either instantiated
    directly from a PDB file (freesasa_structure_from_pdb()) or atom
    by atom (freesasa_structure_add_atom()).

    @defgroup classifier Classifier

    Interface for classifying atoms as polar/apolar and determining
    their radius based on atom name and residue name.

    @defgroup selection Selection

    Interface for selecting a group of atoms and integrating their area.

    @defgroup deprecated Deprecated

    Legacy functions and datatypes from FreeSASA 1.x. Kept because
    they still work although they have been replaced by other
    functions. Can disappear at any time in the future.
 */

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @brief The FreeSASA algorithms. @ingroup core */
enum freesasa_algorithm {
    FREESASA_LEE_RICHARDS, /**< Lee & Richards' algorithm. */
    FREESASA_SHRAKE_RUPLEY /**< Shrake & Rupley's algorithm. */
};

#ifndef __cplusplus
typedef enum freesasa_algorithm freesasa_algorithm;
#endif

/**
   @brief Verbosity levels.
   @see freesasa_set_verbosity()
   @see freesasa_get_verbosity()
 */
enum freesasa_verbosity {
    FREESASA_V_NORMAL,     /**< Print all errors and warnings. */
    FREESASA_V_NOWARNINGS, /**< Print only errors. */
    FREESASA_V_SILENT,     /**< Print no errors and warnings. */
    FREESASA_V_DEBUG,      /**< Print all errors, warnings and debug messages. */
};

#ifndef __cplusplus
typedef enum freesasa_verbosity freesasa_verbosity;
#endif

/* Default parameters */
#define FREESASA_DEF_ALGORITHM FREESASA_LEE_RICHARDS /**< Default algorithm @ingroup core. */
#define FREESASA_DEF_PROBE_RADIUS 1.4                /**< Default probe radius (in Ångström) @ingroup core. */
#define FREESASA_DEF_SR_N 100                        /**< Default number of test points in S&R @ingroup core. */
#define FREESASA_DEF_LR_N 20                         /**< Default number of slices per atom in L&R @ingroup core. */

/**
   @brief Default ::freesasa_classifier
   @ingroup core
 */
#define freesasa_default_classifier freesasa_protor_classifier

/**
   String returned by freesasa_structure_classifier_name() when
   structure was initialized using several different
   classifiers.

   @ingroup structure
 */
#define FREESASA_CONFLICTING_CLASSIFIERS "conflicting-classifiers"

/**
   @brief Default number of threads.

   Value will depend on if library was compiled with or without thread support.
   (2 with threads, 1 without).

   @ingroup core
 */
const extern int FREESASA_DEF_NUMBER_THREADS;

/**
   @brief Error codes.

   Can rely upon FREESASA_SUCCESS being 0 and the errors
   having negative numbers.
 */
enum freesasa_error_codes {
    FREESASA_SUCCESS = 0, /**< All is ok (value will always be zero). */
    FREESASA_FAIL = -1,   /**< Something went seriously wrong (value will always be negative). */
    FREESASA_WARN = -2,   /**< Something went wrong, but results might still be meaningful (value will always be negative). */
};

/**
   @brief Atom classes

   Atoms can be of the classes Apolar, Polar or Unknown.
   @ingroup classifier
 */
enum freesasa_atom_class {
    FREESASA_ATOM_APOLAR = 0,
    FREESASA_ATOM_POLAR = 1,
    FREESASA_ATOM_UNKNOWN = 2,
};

#ifndef __cplusplus
typedef enum freesasa_atom_class freesasa_atom_class;
#endif

/**
   @brief Options for reading structure from PDB

   To be combined in options bitfield in freesasa_structure_from_pdb(),
   freesasa_structure_array() and freesasa_structure_add_atom_wopt().
   See documentation for each function for which options are applicable.

   @ingroup structure
 */
enum freesasa_structure_options {
    FREESASA_INCLUDE_HETATM = 1,             /**< Include HETATM entries. */
    FREESASA_INCLUDE_HYDROGEN = 1 << 2,      /**< Include hydrogen atoms. */
    FREESASA_SEPARATE_MODELS = 1 << 3,       /**< Read MODELs as separate structures. */
    FREESASA_SEPARATE_CHAINS = 1 << 4,       /**< Read separate chains as separate structures. */
    FREESASA_JOIN_MODELS = 1 << 5,           /**< Read MODELs as part of one big structure. */
    FREESASA_HALT_AT_UNKNOWN = 1 << 6,       /**< Halt reading when unknown atom is encountered. */
    FREESASA_SKIP_UNKNOWN = 1 << 7,          /**< Skip atom when unknown atom is encountered. */
    FREESASA_RADIUS_FROM_OCCUPANCY = 1 << 8, /**< Read atom radius from occupancy field. */
};

/**
   @brief Output options

   Controls output format, can be combined in options bitfield in freesasa_tree_export()

   @ingroup node
 */
enum freesasa_output_options {
    FREESASA_OUTPUT_ATOM = 1,           /**< Output data for atoms, residues, chains and structure. */
    FREESASA_OUTPUT_RESIDUE = 1 << 2,   /**< Output data for residues, chains and structure. */
    FREESASA_OUTPUT_CHAIN = 1 << 3,     /**< Output data for chains and structure. */
    FREESASA_OUTPUT_STRUCTURE = 1 << 4, /**< Output data only for the whole structure. */
    FREESASA_LOG = 1 << 5,              /**< Simple plain text results. */
    FREESASA_RSA = 1 << 6,              /**< RSA output (not affected by atom, residue, etc above). */
    FREESASA_JSON = 1 << 7,             /**< JSON output. */
    FREESASA_XML = 1 << 8,              /**< XML output. */
    FREESASA_PDB = 1 << 9,              /**< PDB output (with B-factors replaced by SASA values, and occupancy by radius). */
    FREESASA_RES = 1 << 10,             /**< A list of the integrated SASA of each residue type. */
    FREESASA_SEQ = 1 << 11,             /**< The SASA of each residue in the sequence. */
    FREESASA_CIF = 1 << 12,             /**< CIF output with SASA values and SASA radius appended */

    /**
       Don't output relative areas, for example if structure has
       manually set radii, invalidating reference values
     */
    FREESASA_OUTPUT_SKIP_REL = 1 << 13,
};

/**
   The maximum length of a selection name
   @see freesasa_select_area()
   @ingroup selection
 */
#define FREESASA_MAX_SELECTION_NAME 50

/**
   Struct to store parameters for SASA calculation
   @ingroup core
 */
struct freesasa_parameters {
    freesasa_algorithm alg;     /**< Algorithm. */
    double probe_radius;        /**< Probe radius (in Ångström). */
    int shrake_rupley_n_points; /**< Number of test points in S&R calculation. */
    int lee_richards_n_slices;  /**< Number of slices per atom in L&R calculation. */
    int n_threads;              /**< Number of threads to use, if compiled with thread-support. */
};

#ifndef __cplusplus
typedef struct freesasa_parameters freesasa_parameters;
#endif

/**
   The default parameters for FreeSASA.
   @ingroup core
 */
extern const freesasa_parameters freesasa_default_parameters;

/**
   @brief Struct for structure object.

   The struct includes the coordinates and radius of each atom, and
   its name, residue-name, etc. If it was initiated from a PDB file
   enough info will be stored so that a PDB-file can be printed using
   the original one as template.

   @ingroup structure
 */
typedef struct freesasa_structure freesasa_structure;

/**
   Struct to store results of SASA calculation

   @ingroup core
 */
struct freesasa_result {
    double total;                   /**< Total SASA in Ångström^2. */
    double *sasa;                   /**< SASA of each atom in Ångström^2. */
    int n_atoms;                    /**< Number of atoms. */
    freesasa_parameters parameters; /**< Parameters used when generating result. */
};

#ifndef __cplusplus
typedef struct freesasa_result freesasa_result;
#endif

/**
   Struct to store integrated SASA values for either a full structure
   or a subset thereof.

   Use freesasa_result_classes() to turn a
   ::freesasa_result into a ::freesasa_nodearea. Each
   ::freesasa_node is associated with a
   ::freesasa_nodearea.

   @ingroup node
 */
struct freesasa_nodearea {
    const char *name;  /**< Name of substructure. */
    double total;      /**< Total SASA. */
    double main_chain; /**< Main-chain/Backbone SASA. */
    double side_chain; /**< Side-chain SASA. */
    double polar;      /**< Polar SASA. */
    double apolar;     /**< Apolar SASA. */
    double unknown;    /**< SASA of unknown class (neither polar nor apolar). */
};

#ifndef __cplusplus
typedef struct freesasa_nodearea freesasa_nodearea;
#endif

/**
   @brief Node types
   @ingroup node
 */
enum freesasa_nodetype {
    FREESASA_NODE_ATOM,      /**< Atom node. */
    FREESASA_NODE_RESIDUE,   /**< Residue node. */
    FREESASA_NODE_CHAIN,     /**< Chain node. */
    FREESASA_NODE_STRUCTURE, /**< Structure node. */
    FREESASA_NODE_RESULT,    /**< Result node, wraps results for one or more related structures. */
    FREESASA_NODE_ROOT,      /**< Root node, wraps one or more unrelated results. */
    FREESASA_NODE_NONE       /**< for specifying not a valid node. */
};

#ifndef __cplusplus
typedef enum freesasa_nodetype freesasa_nodetype;
#endif

/**
   Struct to store data about a mmCIF atom site.

   @ingroup structure
 */
struct freesasa_cif_atom {
    const char *group_PDB;
    const char auth_asym_id;
    const char *auth_seq_id;
    const char *pdbx_PDB_ins_code;
    const char *auth_comp_id;
    const char *auth_atom_id;
    const char *label_alt_id;
    const char *type_symbol;
    const double Cartn_x;
    const double Cartn_y;
    const double Cartn_z;
};

#ifndef __cplusplus
typedef struct freesasa_cif_atom freesasa_cif_atom;
#endif

/**
   @brief Result node

   A node representing calculation results for a structure, chain,
   residue or atom in a structure (see @ref node).

   @ingroup node
 */
typedef struct freesasa_node freesasa_node;

/**
   @brief Selection struct

   Struct to store a selection generated by freesasa_selection_new().

   @ingroup selection
 */
typedef struct freesasa_selection freesasa_selection;

/**
   @brief Classifier struct

   Struct that can be used to determine classes (polar/apolar) and
   radii of atoms. Initiated from freesasa_classifier_from_file().
   The classifiers ::freesasa_default_classifier, ::freesasa_protor_classifier,
   ::freesasa_naccess_classifier and ::freesasa_classifier are const
   classifiers that can be used directly.

   @ingroup classifier
 */
typedef struct freesasa_classifier freesasa_classifier;

/**
   @brief ProtOr classifier.

   Classifier using ProtOr radii and classes.

   @ingroup classifier
 */
extern const freesasa_classifier freesasa_protor_classifier;

/**
   @brief NACCESS classifier.

   Classifier using NACCESS radii and classes.

   @ingroup classifier
 */
extern const freesasa_classifier freesasa_naccess_classifier;

/**
   @brief OONS classifier.

   Classifier using OONS radii and classes.

   @ingroup classifier
 */
extern const freesasa_classifier freesasa_oons_classifier;

/**
    Calculates SASA based on a given structure.

    This function allows direct access to the results array,
    for most users freesasa_calc_tree() is more appropriate.

    Return value is dynamically allocated, should be freed with
    freesasa_result_free().

    @param structure The structure
    @param parameters Parameters for the calculation, if `NULL`
      defaults are used.

    @return The result of the calculation, `NULL` if something went wrong.

    @ingroup core
 */
freesasa_result *
freesasa_calc_structure(const freesasa_structure *structure,
                        const freesasa_parameters *parameters);

/**
    Calculates SASA based on a given set of coordinates and radii.

    Return value is dynamically allocated, should be freed with
    freesasa_result_free().

    @param xyz Array of coordinates in the form x1,y1,z1,x2,y2,z2,...,xn,yn,zn.
    @param radii Radii, this array should have n elements..
    @param n Number of coordinates (i.e. xyz has size 3*n, radii size n).
    @param parameters Parameters for the calculation, if `NULL`
      defaults are used.

    @return The result of the calculation, `NULL` if something went wrong.

    @ingroup core
 */
freesasa_result *
freesasa_calc_coord(const double *xyz,
                    const double *radii,
                    int n,
                    const freesasa_parameters *parameters);

/**
    Calculates SASA for a structure and returns as a tree of
    ::freesasa_node.

    Return value is dynamically allocated, should be freed with
    freesasa_node_free()

    @param structure A structure
    @param parameters Parameters for the calculation, if `NULL`,
      defaults are used.
    @param name Name of input structure to be used in output.

    @return The result of the calculation, `NULL` if something went wrong.

    @ingroup core
 */
freesasa_node *
freesasa_calc_tree(const freesasa_structure *structure,
                   const freesasa_parameters *parameters,
                   const char *name);

/**
    Results by classes.

    Adds up the SASA of Polar/Apolar/Unknown atoms, and
    main-chain/side-chain atoms for the whole protein. Uses the
    classes defined by the classifier used when generating the
    structure.

    @param structure The structure the results are based on
    @param result The results
    @return A struct with all the results.

    @ingroup core
 */
freesasa_nodearea
freesasa_result_classes(const freesasa_structure *structure,
                        const freesasa_result *result);

/**
    Frees a ::freesasa_result object.

    @param result the object to be freed.

    @ingroup core
 */
void freesasa_result_free(freesasa_result *result);

/**
    Generate a classifier from a config-file.

    Input file format described in @ref Config-file

    Return value is dynamically allocated, should be freed with
    freesasa_classifier_free().

    @param file File containing configuration
    @return The generated classifier. `NULL` if there were problems
      parsing or reading the file or memory allocation problem.

    @see @ref Config-file

    @ingroup classifier
 */
freesasa_classifier *
freesasa_classifier_from_file(FILE *file);

/**
    Frees a classifier object

    @param classifier The classifier.

    @ingroup classifier
 */
void freesasa_classifier_free(freesasa_classifier *classifier);

/**
    Use a classifier to determine the radius of a given atom.

    @param classifier The classifier.
    @param res_name The residue name (ALA/VAL/U/C/...)
    @param atom_name The atom name (CA/N/CB/...)
    @return The radius, negative if atom unknown.

    @ingroup classifier
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

    @ingroup classifier
 */
freesasa_atom_class
freesasa_classifier_class(const freesasa_classifier *classifier,
                          const char *res_name,
                          const char *atom_name);

/**
    Names for ::freesasa_atom_class.

    @param atom_class The class.
    @return Name of class.

    @ingroup classifier
 */
const char *
freesasa_classifier_class2str(freesasa_atom_class atom_class);

/**
    The name of a classifier.

    @param classifier The classifier.
    @return The name of the classifier.

    @ingroup classifier
 */
const char *
freesasa_classifier_name(const freesasa_classifier *classifier);

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

    The return value should be freed with freesasa_selection_free().

    @see @ref Selection

    @param command The selection
    @param structure The structure to select from
    @param result The results to integrate
    @return The selection. `NULL` if something went wrong. Use
      freesasa_selection_name(), freesasa_selection_command(),
      freesasa_selection_area() and freesasa_selection_n_atoms() to
      access results of selection.

    @ingroup selection
*/
freesasa_selection *
freesasa_selection_new(const char *command,
                       const freesasa_structure *structure,
                       const freesasa_result *result);

/**
    Free selection.

    @param selection The selection

    @ingroup selection
 */
void freesasa_selection_free(freesasa_selection *selection);

/**
    Name of the selection

    @param selection The selection
    @return the name

    @ingroup selection
 */
const char *
freesasa_selection_name(const freesasa_selection *selection);

/**
    Command that was used to generate the selection

    @param selection The selection
    @return The command

    @ingroup selection
 */
const char *
freesasa_selection_command(const freesasa_selection *selection);

/**
    Area of the selection

    @param selection The selection
    @return The area

    @ingroup selection
 */
double
freesasa_selection_area(const freesasa_selection *selection);

/**
    Number of atoms that matched the selection

    @param selection The selection
    @return Number of atoms

    @ingroup selection
 */
int freesasa_selection_n_atoms(const freesasa_selection *selection);

/**
    Set the global verbosity level.

    @param v the verbosity level
    @return ::FREESASA_SUCCESS. If v is invalid ::FREESASA_FAIL.
    @see freesasa_verbosity

    @ingroup core
 */
int freesasa_set_verbosity(freesasa_verbosity v);

/**
    Get the current verbosity level

    @return the verbosity level.

    @ingroup core
 */
freesasa_verbosity
freesasa_get_verbosity(void);

/**
    Set where to write errors.

    By default `stderr` is used, this function can be called to redirect
    error output elsewhere.

    @param err The file to write to. If `NULL`, `stderr` will be used.

    @ingroup core
 */
void freesasa_set_err_out(FILE *err);

/**
    Get pointer to error file.

    `NULL` means `stderr` is used.

    @return The error file.

    @ingroup core
 */
FILE *
freesasa_get_err_out(void);

/**
    Allocate empty structure.

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @return Empty structure. `NULL` if memory allocation failure.

    @ingroup structure
 */
freesasa_structure *
freesasa_structure_new(void);

/**
    Free structure.

    @param structure The structure to free.

    @ingroup structure
 */
void freesasa_structure_free(freesasa_structure *structure);

/**
    Init structure with coordinates from pdb-file.

    Reads in a PDB-file and generates a structure object.
    Automatically skips hydrogens and HETATM. If an atom has
    alternative coordinates, only the first alternative is used. If a
    file has more than one MODEL (as in NMR structures) only the
    first model is used. The provided classifier is used to determine
    the radius of each atom, if the atom is not recognized the element
    of the atom is guessed, and that element's VdW radius used. If
    this fails its radius is set 0, which means that it won't
    contribute to SASA, but a radius from another source can be set
    through freesasa_structure_set_radius(). All these behaviors can
    be modified through the `options` bitfield argument:

      - 0: Default behavior

      - ::FREESASA_INCLUDE_HYDROGEN: Include hydrogen atoms.

      - ::FREESASA_INCLUDE_HETATM: Include HETATM.

      - ::FREESASA_JOIN_MODELS: Join models.

      - ::FREESASA_SKIP_UNKNOWN: Skip unknown atoms.

      - ::FREESASA_HALT_AT_UNKNOWN: Halt at unknown atom and return
        `NULL`. Overrides ::FREESASA_SKIP_UNKNOWN.

      - ::FREESASA_RADIUS_FROM_OCCUPANCY: Read atomic radii from
         Occupancy field in PDB file.

    If a more fine-grained control over which atoms to include is
    needed, the PDB-file needs to be modified before calling this
    function, or atoms can be added manually one by one using
    freesasa_structure_add_atom() or
    freesasa_structure_add_atom_wopt().

    Return value is dynamically allocated, should be freed with
    freesasa_structure_free().

    @param pdb A PDB file

    @param classifier A freesasa_classifier to determine radius of
      atom. If `NULL` default classifier is used.

    @param options A bitfield to determine what atoms to include and what to do
      when atoms are not recognized by classifier.

    @return The generated structure. Returns `NULL` and prints error
      if input is invalid or  memory allocation failure.

    @ingroup structure
 */
freesasa_structure *
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
      chain, respectively (will return `NULL` if neither is
      specified). See freesasa_structure_from_pdb() for documentation
      on options for deciding what atoms to include
      (::FREESASA_JOIN_MODELS is not supported here).
    @return Array of structures. Prints error message(s) and returns
      `NULL` if there were problems reading input, if invalid value of
      `options`, or upon a memory allocation failure.

    @ingroup structure
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

    @ingroup structure
 */
int freesasa_structure_add_atom(freesasa_structure *structure,
                                const char *atom_name,
                                const char *residue_name,
                                const char *residue_number,
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

      - `options & ::FREESASA_SKIP_UNKNOWN == 1` skip unknown atoms,
        return ::FREESASA_WARN.

      - `options & ::FREESASA_HALT_AT_UNKNOWN == 1` skip unknown atoms,
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
    @param classifier A ::freesasa_classifier to determine radius of atom and to
      decide if to keep atom or not (see options).
    @param options A bitfield to determine what to do with unknown atoms (see above).

    @return ::FREESASA_SUCCESS on normal execution. ::FREESASA_FAIL if
       if memory allocation fails or if halting at unknown
       atom. ::FREESASA_WARN if skipping atom.

    @ingroup structure
 */
int freesasa_structure_add_atom_wopt(freesasa_structure *structure,
                                     const char *atom_name,
                                     const char *residue_name,
                                     const char *residue_number,
                                     char chain_label,
                                     double x, double y, double z,
                                     const freesasa_classifier *classifier,
                                     int options);
/**
    Add atoms from a mmCIF file to a structure

    @param structure The structure to add to.
    @param atom An atom site from a mmCIF file
    @param classifier A ::freesasa_classifier to determine radius of atom and to
      decide if to keep atom or not (see options).
    @param options Structure options as in freesasa_structure_add_atom_wopt()

    @return ::FREESASA_SUCCESS on normal execution. ::FREESASA_FAIL if
       if memory allocation fails or if halting at unknown
       atom. ::FREESASA_WARN if skipping atom.

    @ingroup structure
 */
int freesasa_structure_add_cif_atom(freesasa_structure *structure,
                                    freesasa_cif_atom *atom,
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
    @param chains String of chain labels (e.g. `"AB"`)
    @param classifier A classifier to use to build the new structure
    @param options Structure options as in freesasa_structure_add_atom_wopt()

    @return A new structure consisting only of the specified
    chains. Returns `NULL` if one or more of the requested chains don't
    match any in the input structure or if memory allocation fails.

    @ingroup structure
 */
freesasa_structure *
freesasa_structure_get_chains(const freesasa_structure *structure,
                              const char *chains,
                              const freesasa_classifier *classifier,
                              int options);

/**
    Get string listing all chains in structure.

    @param structure The structure.
    @return String with all chain labels in structure (`"A"`, `"ABC"`, etc).

    @ingroup structure
 */
const char *
freesasa_structure_chain_labels(const freesasa_structure *structure);

/**
    Get number of atoms.

    @param structure The structure.
    @return Number of atoms.

    @ingroup structure
 */
int freesasa_structure_n(const freesasa_structure *structure);

/**
    Get number of residues.

    Calculated crudely by determining the number of unique
    combinations of residue number and chain label contained in the
    structure. If residues are mingled i.e. atoms of the same residue
    are in non-contiguous regions of the file, this function might be
    off.

    @param structure A structure.
    @return Number of residues.

    @ingroup structure
 */
int freesasa_structure_n_residues(const freesasa_structure *structure);

/**
    Get number of chains.

    @param structure A structure.
    @return The number of chains in the structure.

    @ingroup structure
 */
int freesasa_structure_n_chains(const freesasa_structure *structure);

/**
    Returns a pointer to an array of the radii of each atom.

    @param structure The structure.
    @return Array of radii. If `NULL` structure has not been properly
      initialized.

    @ingroup structure
 */
const double *
freesasa_structure_radius(const freesasa_structure *structure);

/**
    Override the radii set when the structure was initialized.

    Makes a copy of the provided array.

    @param structure The structure.
    @param radii An array of radii, should have same dimension
      as the number of atoms in the structure.

    @ingroup structure
 */
void freesasa_structure_set_radius(freesasa_structure *structure,
                                   const double *radii);

/**
    Get atom name.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.

    @return Atom name in the form `" CA "`, `" OXT"`, etc, if
       structure was initialized from a PDB file, or in whatever form
       it was added through freesasa_structure_add_atom() or
       freesasa_structure_add_atom_wopt().

    @ingroup structure
 */
const char *
freesasa_structure_atom_name(const freesasa_structure *structure,
                             int i);

/**
    Get residue name.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Residue name in the form `"ALA"`, `"PHE"`, `"  C"`, etc, if
       structure was initialized from a PDB file, or in whatever form
       it was added through freesasa_structure_add_atom() or
       freesasa_structure_add_atom_wopt().

    @ingroup structure
 */
const char *
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

    @ingroup structure
 */
const char *
freesasa_structure_atom_res_number(const freesasa_structure *structure,
                                   int i);

/**
    Get chain label.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Chain label (`'A'`, `'B'`, etc.)

    @ingroup structure
 */
char freesasa_structure_atom_chain(const freesasa_structure *structure,
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

    @ingroup structure
 */
const char *
freesasa_structure_atom_symbol(const freesasa_structure *structure,
                               int i);

/**
    Get atom radius.

    Asserts that index i is within bounds.

    @param structure The structure.
    @param i Atom index.
    @return Atom radius.

    @ingroup structure
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

    @ingroup structure
 */
void freesasa_structure_atom_set_radius(freesasa_structure *structure,
                                        int i,
                                        double radius);

/**
    Get name of residue.

    @param structure The structure.
    @param r_i Residue index (in whole structure)
    @return Name of residue

    @ingroup structure
 */
const char *
freesasa_structure_residue_name(const freesasa_structure *structure,
                                int r_i);

/**
    Get residue number.

    @param structure The structure.
    @param r_i Residue index (in whole structure).
    @return Residue number as string.

    @ingroup structure
 */
const char *
freesasa_structure_residue_number(const freesasa_structure *structure,
                                  int r_i);
/**
    Get chain residue belongs to.

    @param structure The structure.
    @param r_i Residue index (in whole structure).
    @return Chain label.

    @ingroup structure
 */
char freesasa_structure_residue_chain(const freesasa_structure *structure,
                                      int r_i);

/**
    Get model number for structure.

    Useful if structure was generated with freesasa_structure_array().

    @param structure The structure.
    @return The model number. Will be 1 if MODEL not specified in PDB
      input.

    @ingroup structure
 */
int freesasa_structure_model(const freesasa_structure *structure);

/**
    Get array of coordinates.

    Size of array is 3*N, order of coordinates `x1, y1, z1, ...`.

    @param structure The structure.
    @return Array of coordinates. `NULL` if structure empty. Size can be
      accessed through freesasa_structure_n() (multiply by three).

    @ingroup structure
 */
const double *
freesasa_structure_coord_array(const freesasa_structure *structure);

/**
    Get indices of first and last atoms of a residue

    @param structure A structure.
    @param r_i Residue index.
    @param first First atom of residue `r_i` will be stored here.
    @param last Last atom of residue `r_i` will be stored here.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if index `r_i` is invalid.

    @ingroup structure
 */
int freesasa_structure_residue_atoms(const freesasa_structure *structure,
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

    @ingroup structure
 */
int freesasa_structure_chain_atoms(const freesasa_structure *structure,
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

    @ingroup structure
 */
int freesasa_structure_chain_residues(const freesasa_structure *structure,
                                      char chain,
                                      int *first,
                                      int *last);

/**
    Name of classifier used to generate structure.

    @param structure A structure.
    @return Name of classifier. Name will equal
      ::FREESASA_CONFLICTING_CLASSIFIERS if several different
      classifiers were used.

    @ingroup structure
 */
const char *
freesasa_structure_classifier_name(const freesasa_structure *structure);

/**
    Generates empty ::freesasa_node of type ::FREESASA_NODE_ROOT.

    To be populated by freesasa_tree_add_result().

    The return value is dynamically allocated and should be freed
    using freesasa_node_free().

    @return A ::freesasa_node. `NULL` if memory allocation fails.

    @ingroup node
 */
freesasa_node *
freesasa_tree_new(void);

/**
    Init tree based on result and structure.

    The return value is dynamically allocated and should be freed
    using freesasa_node_free().

    @param result A result.
    @param structure A structure.
    @param name Name of the results (typically filename from
      which structure is derived)

    @return The root node of the tree. `NULL` if memory allocation
      fails.

    @ingroup node
*/
freesasa_node *
freesasa_tree_init(const freesasa_result *result,
                   const freesasa_structure *structure,
                   const char *name);

/**
    Add a new set of results to a tree.

    Tree should first be initiated with freesasa_calc_tree(),
    freesasa_tree_new() or freesasa_tree_init().

    @param tree Node of type ::FREESASA_NODE_ROOT. Tree to add results
      to.
    @param result SASA values for the structure
    @param structure The structure the results are based on
    @param name The name to use for the result
    @return ::FREESASA_SUCCESS upon success. ::FREESASA_FAIL if memory
      allocation fails.

    @ingroup node
 */
int freesasa_tree_add_result(freesasa_node *tree,
                             const freesasa_result *result,
                             const freesasa_structure *structure,
                             const char *name);

/**
    Join two trees.

    Allows joining several calculations into one output file.

    @param tree1 Node of type ::FREESASA_NODE_ROOT. The joint tree
      will be stored here.
    @param tree2 Node of type ::FREESASA_NODE_ROOT. Will be added to
      tree1, and then changed to `NULL`, since ownership of its contents
      have been transferred to tree1.
    @return ::FREESASA_SUCCESS.

    @ingroup node
 */
int freesasa_tree_join(freesasa_node *tree1,
                       freesasa_node **tree2);

/**
    Outputs result in format specified by options.

    @param output Output file.
    @param root Structure tree containing results. Node of type ::FREESASA_NODE_ROOT.
    @param options Bitfield specifying output format, see
      ::freesasa_output_options.
    @return ::FREESASA_SUCCESS upon success. ::FREESASA_FAIL if there
      was an error (see messages).

    @ingroup node
*/
int freesasa_tree_export(FILE *output,
                         freesasa_node *root,
                         int options);

/**
    Free tree.

    Will not free anything if the node is not a root node.

    @param root Node of type ::FREESASA_NODE_ROOT
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if the node has a
      parent.

    @ingroup node
 */
int freesasa_node_free(freesasa_node *root);

/**
    The ::freesasa_nodearea of all atoms belonging to a node.

    @param node The node.
    @return The area. `NULL` if no area has been attached to this node.

    @ingroup node
 */
const freesasa_nodearea *
freesasa_node_area(const freesasa_node *node);

/**
    The children of a node.

    Use freesasa_node_next() to access next sibling.

    @param node The node.
    @return Pointer to the first child of a node. `NULL` if the node has no
      children.

    @ingroup node
 */
freesasa_node *
freesasa_node_children(freesasa_node *node);

/**
    Next sibling of a node.

    @param node The node.
    @return The next node, `NULL` if this is the last node.

    @ingroup node
 */
freesasa_node *
freesasa_node_next(freesasa_node *node);

/**
    The parent of a node.

    @param node The node.
    @return The parent node. `NULL` if the node has no parent.

    @ingroup node
 */
freesasa_node *
freesasa_node_parent(freesasa_node *node);

/**
    The type of a node.

    @param node The node.
    @return The type.

    @ingroup node
 */
freesasa_nodetype
freesasa_node_type(const freesasa_node *node);

/**
    The name of a node.

    The node types will have the following names:
    - Atom: atom name, i.e. `" CA "`, `" OXT"`, etc.
    - Residue: residue name, i.e. `"ALA"`, `"ARG"`, etc.
    - Chain: chain label, i.e. `"A"`, `"B"`, etc.
    - Structure: string of all chain labels in the molecule, i.e. `"A"`, `"ABC"`, etc
    - Result: name of input (most often input filename or `"stdin"`)
    - Root: `NULL`

    @param node The node.
    @return The name. `NULL` if the node has no name.

    @ingroup node
 */
const char *
freesasa_node_name(const freesasa_node *node);

/**
    The name of the classifier used to generate the node.

    @param node A node of type ::FREESASA_NODE_RESULT.
    @return The name of the classifier

    @ingroup node
 */
const char *
freesasa_node_classified_by(const freesasa_node *node);

/**
    Is atom polar.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return 1 if polar, 0 else.

    @ingroup node
 */
int freesasa_node_atom_is_polar(const freesasa_node *node);

/**
    Does atom belong to the main chain/backbone.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return 1 if mainchain, 0 else.

    @ingroup node
 */
int freesasa_node_atom_is_mainchain(const freesasa_node *node);

/**
    Atom radius.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return The radius.

    @ingroup node
 */
double
freesasa_node_atom_radius(const freesasa_node *node);

/**
    Line in PDB atom was generated from.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return The line. `NULL` if atom wasn't taken from PDB file.

    @ingroup node
 */
const char *
freesasa_node_atom_pdb_line(const freesasa_node *node);

/**
    Atom residue number.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return The residue sequence number that this atom is a part of.

    @ingroup node
 */
const char *
freesasa_node_atom_residue_number(const freesasa_node *node);

/**
    Atom residue name.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return The residue 3-char name this atom is a part of.

    @ingroup node
 */
const char *
freesasa_node_atom_residue_name(const freesasa_node *node);

/**
    Atom chain.

    @param node A node of type ::FREESASA_NODE_ATOM.
    @return The chain this atom is a part of.

    @ingroup node
 */
char freesasa_node_atom_chain(const freesasa_node *node);

/**
    Residue number.

    @param node A node of type ::FREESASA_NODE_RESIDUE.
    @return String with residue number.

    @ingroup node
 */
const char *
freesasa_node_residue_number(const freesasa_node *node);

/**
    Number of atoms in a residue.

    @param node A node of type ::FREESASA_NODE_RESIDUE.
    @return Number of atoms.

    @ingroup node
 */
int freesasa_node_residue_n_atoms(const freesasa_node *node);

/**
    The reference area for a node from the classifier used to
    generate the tree.

    @param node A node of type ::FREESASA_NODE_RESIDUE.
    @return The reference area. `NULL` if area not available.

    @ingroup node
 */
const freesasa_nodearea *
freesasa_node_residue_reference(const freesasa_node *node);

/**
    The number of residues in a chain.

    @param node A node of type ::FREESASA_NODE_CHAIN.
    @return Number of residues.

    @ingroup node
 */
int freesasa_node_chain_n_residues(const freesasa_node *node);

/**
    The number of chains in a structure.

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @return Number of chains.

    @ingroup node
 */
int freesasa_node_structure_n_chains(const freesasa_node *node);

/**
    The number of atoms in a structure.

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @return Number of atoms.

    @ingroup node
 */
int freesasa_node_structure_n_atoms(const freesasa_node *node);

/**
    All chain labels in a structure.

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @return Chain labels as null-terminated string.

    @ingroup node
 */
const char *
freesasa_node_structure_chain_labels(const freesasa_node *node);

/**
    Model number of a structure (from input PDB)

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @return Model number.

    @ingroup node
 */
int freesasa_node_structure_model(const freesasa_node *node);

/**
    Raw results for a structure

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @return The results.

    @ingroup node
 */
const freesasa_result *
freesasa_node_structure_result(const freesasa_node *node);

/**
    Selection results for a structure

    Generated using freesasa_node_structure_add_selection().

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @return A null-terminated array of pointers to selections. `NULL` if
      no selections were associated with structure.

    @ingroup node
 */
const freesasa_selection **
freesasa_node_structure_selections(const freesasa_node *node);

/**
    Add a selection result to a structure node

    The selection is cloned, so the user can call
    freesasa_selection_free() on the provided selection at the time of
    their chosing.

    @param node A node of type ::FREESASA_NODE_STRUCTURE.
    @param selection A selection.
    @return ::FREESASA_SUCCESS. ::FREESASA_FAIL if cloning fails
      (i.e. memory allocation failure).

    @ingroup node
 */
int freesasa_node_structure_add_selection(freesasa_node *node,
                                          const freesasa_selection *selection);
/**
    Parameter values used to calculate result.

    @param node A node of type ::FREESASA_NODE_RESULT
    @return The parameters.

    @ingroup node
 */
const freesasa_parameters *
freesasa_node_result_parameters(const freesasa_node *node);

/* Deprecated functions below, from 1.x API */

/**
    Get area of a selection.

    @deprecated Use freesasa_selection_new() instead.

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

    @ingroup deprecated
 */
int freesasa_select_area(const char *command,
                         char *name,
                         double *area,
                         const freesasa_structure *structure,
                         const freesasa_result *result);

void freesasa_structure_set_model(freesasa_structure *structure,
                                  int model);

#ifdef __cplusplus
}
#endif

#endif /* FREESASA_H */
