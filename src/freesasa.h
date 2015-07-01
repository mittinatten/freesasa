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

/**
   Type used to store parameters and results of FreeSASA
   calculations. 
 */
typedef struct freesasa freesasa;

/**
   Type used to store SASA values for a groups of atoms and strings
   describing those groups of atoms. The parameter n specifies how
   many elements the arrays have.
 */
typedef struct {
    double *value;
    char **string;
    int n;
} freesasa_strvp;

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
                         
/////////////////////////////////////////
// Initialization, freeing and copying //
/////////////////////////////////////////

#ifdef __cplusplus
extern "C"{
#endif

/**
    Initiate freesasa-object

    Allocates empty freesasa object with default parameters. Must be
    called before calculation. 

    @return The created object
 */
    freesasa* freesasa_new(void);

/**
    Frees resources allocated to a freesasa object. 

    @param s a ::freesasa-object
 */
    void freesasa_free(freesasa *s);

/**
    Copy parameters
    
    Copies all algorithm parameters. Data from calculations are not
    copied, i.e. this can be used for setting up several calculations
    with identical parameters for different proteins.
    
    @param target The target configuration
    @param source The source configuration
 */
    void freesasa_copy_param(freesasa *target, const freesasa *source);


//////////////////////////
// Perform calculations //
//////////////////////////

/**
    Calculate SASA from coordinates and radii. 

    Can be run on any set of spheres. Coordinates and radii will not
    be stored. This will only generate a total SASA, not polar or
    apolar, since only radii are provided, not atom types.

    @see freesasa_calc_pdb()
    @see freesasa_calc_atoms()
    @param s a ::freesasa-object
    @param coord An array of coordinates for sphere centers
      (x1,y1,z1,x2,y2,z2,...)
    @param r An array of radii for the spheres
    @param n number of spheres
    @return ::FREESASA_SUCCESS upon successful calculation, prints and
    error and returns ::FREESASA_FAIL else. 
 */
    int freesasa_calc_coord(freesasa  *s, const double *coord,
                            const double *r, int n);

/**
    Calculate SASA from PDB-file.
    
    Reads PDB-file and calculates SASA. HETATM records are
    ignored. Results stored in parameter s. If s is not initialized
    default values are used, these are stored in s.  If the object has
    been used in calculations previously, the results from these will
    be over-written. 

    @see freesasa_calc_atoms()
    @see freesasa_calc_coord()
    @param s a ::freesasa-object
    @param pdb_file PDB-file for input
    @return ::FREESASA_SUCCESS if calculation successful, prints an
    error and returns ::FREESASA_FAIL if not.
 */
    int freesasa_calc_pdb(freesasa *s, FILE *pdb_file);

/**
    Calculate SASA from a set of atoms.
    
    Similar to freesasa_calc_pdb(), but takes array of coordinates,
    and residue and atom names separately. The format for resnames and
    atomnames is the same as used in freesasa_radius(). These are
    necessary to determine radii and classify atoms
    (polar/apolar/etc). The same object can be reused for repeated
    calculations on different structures (i.e. the necessary variables
    are reset when it is called). 

    @see freesasa_calc_pdb()
    @see freesasa_calc_coord()
    @param s a ::freesasa-object
    @param coord Array of atom coordinates (x1,y1,z1,x2,y2,z2,...). Should
      have 3*n elements.
    @param resnames Array of strings of the format "ALA", "PHE", etc. 
      Should have n elements.
    @param atomnames Array of strings of the format " CA ", " OXT", etc. 
      Should have n elements.
    @param n number of atoms (>0)

    @return ::FREESASA_SUCCESS if calculation successful.
      ::FREESASA_WARN if atoms can't be classified. ::FREESASA_FAIL if
      if calculations failed.
*/
    int freesasa_calc_atoms(freesasa *s, const double *coord,
                            const char **resnames, 
                            const char **atomnames, int n);

/* Reads pdb-file and calculates radii for each atom. Memory is
   allocated to store them in the array 'r'. The return value is the
   size of the allocated array. Prints error and returns FREESASA_FAIL
   if reading input fails. This can be used if coordinates are linked
   as below and default radii are to be used.  Not properly tested
   yet!
*/
//int freesasa_generate_radii(double **r, FILE *pdb_file);

/**
    Default radius of an atom type. 

    Returns the radius of an atom based on it's type, either according
    to the OONS classification, or the element if OONS class cannot be
    determined. Unknown atom types and hydrogens are assigned radius
    0.0. The residue and atom names are the default names used in PDB
    files. Any whitespace in the standard needs to be included here,
    i.e. "CA" should be " CA ".
        
    @param residue_name Residue name in the PDB format, "ALA", "PHE", etc.
    @param atom_name Atom name, " CA ", " OXT", etc.
    @return Radius of atom in Ångström. 
*/
    double freesasa_radius(const char *residue_name, const char *atom_name);

/**
    Link a set of coordinates to the freesasa object.

    If the linked coordinates are updated freesasa_refresh() can be
    used to recalculate SASA. FreeSASA will not change the
    coordinates. If the freesasa-object has been initalized with a
    PDB file the corresponding memory is released and the results
    erased.
    
    @see freesasa_refresh()
    @param s a ::freesasa-object
    @param coord An array of coordinates for sphere centers (x1,y1,z1,x2,y2,z2,...)
    @param r An array of radii for the spheres
    @param n number of spheres (>0)
    @return At the moment always returns ::FREESASA_SUCCESS
*/
    int freesasa_link_coord(freesasa *s, const double *coord,
                            double *r, int n);

/**
    Recalculate SASA based on linked coordinates.

    Recalculates SASA, based on the assumption that a set of external
    coordinates have been updated elsewhere. 

    @see freesasa_link_coord()
    @param s a ::freesasa-object. Used for parameters, coordinate-link
    and to store results.  
    @return ::FREESASA_SUCCESS upon successful computation. ::FREESASA_FAIL else.
*/
    int freesasa_refresh(freesasa *s);

/**
    The number of atoms in the latest SASA calculation.  
    
    Returns 0 if no coordinates have been linked or no calculation
    has been performed.

    @param s a ::freesasa-object
    @return The number of atoms
*/
    int freesasa_n_atoms(const freesasa *s);

/**
    The number of residues in the protein.

    Requires that a protein, not only raw coordinates, has been
    associated with the provided ::freesasa-object. That is,
    freesasa_calc_pdb() or freesasa_calc_atoms() have been called
    successfully beforehand.

    @param s a ::freesasa-object
    @return The number of residues.
 */
    int freesasa_n_residues(const freesasa *s);

//////////////////////////////
// Settings for calculation //
//////////////////////////////

/**
    Set algorithm. 

    @param s a ::freesasa-object
    @param alg The algorithm to be set.
*/
    void freesasa_set_algorithm(freesasa *s, freesasa_algorithm alg);

/**
    Get algorithm. 

    @param s a ::freesasa-object
    @return The algorithm.
*/
    freesasa_algorithm freesasa_get_algorithm(const freesasa *s);

/**
    Get algorithm as string.
 
    @param s a ::freesasa-object
    @return The algorithm as a string.
*/
    const char* freesasa_algorithm_name(const freesasa *s);

/**
    Set probe radius.

    Sets probe radius for SASA calculations (default
    ::FREESASA_DEF_PROBE_RADIUS = 1.4 Å). If submitted radius is
    invalid, previous value is kept and an error message printed.

    @param s a ::freesasa-object
    @param r Value for probe radius in Ångström.
    @return ::FREESASA_SUCCESS for valid r-values. ::FREESASA_WARN else. 
*/
    int freesasa_set_probe_radius(freesasa *s,double r);

/**
    Get probe radius.

    @param s a ::freesasa-object
    @return Probe radius in Ångström. 
*/
    double freesasa_get_probe_radius(const freesasa *s);
    
/**
    Set number of test points for S&R algorithm.

    Sets number of test points for S&R algorithm (default
    ::FREESASA_DEF_SR_N = 100).
    
    @param s a ::freesasa-object
    @param n Number of points.
    @return ::FREESASA_SUCCESS if n is valid. Else: prints error
    message, returns ::FREESASA_WARN and sets to default value.
*/
    int freesasa_set_sr_points(freesasa *s, int n);

/**
    Get number of test points for S&R algorithm. 

    @param s a ::freesasa-object
    @return Number of points. ::FREESASA_WARN if S&R algorithm not
    selected.
*/
    int freesasa_get_sr_points(const freesasa *s);

/**
    Set L&R slice width.

    Sets slice width d for L&R algorithm in Ångström (default
    FREESASA_DEF_LR_D = 0.25 Å). 

    @param s a ::freesasa-object
    @param d Slice width.

    @return ::FREESASA_SUCCESS if 0 < d <= 2.0. Larger d are accepted
    but the functions returns ::FREESASA_WARN and prints
    warning. Returns ::FREESASA_WARN and sets to default value for d
    <= 0 or when d is not finite.
    
*/
    int freesasa_set_lr_delta(freesasa *s, double d);

/**
    Get L&R slice width.
    
    Get slice width for L&R algorithm in Ångström. 

    @param s a ::freesasa-object
    @return Slice width. Negative value if L&R algorithm not
    selected. 
*/
    double freesasa_get_lr_delta(const freesasa *s);

/**
    Set number of threads.

    Sets the number of threads to use for parallel computations,
    useful for large proteins and high resolution. Minimal gain for
    smaller systems.

    @param s a ::freesasa-object
    @param n Number of threads.
    @return ::FREESASA_SUCCESS if n is valid, uses default value and
    returns ::FREESASA_WARN else.
*/
    int freesasa_set_nthreads(freesasa *s,int n);

/**
    Get number of threads.

    Gets the number of threads used in the calcultion.

    @param s a ::freesasa-object
    @return Number of threads.
*/
    int freesasa_get_nthreads(const freesasa *s);

/**
    Set protein name
    
    Sets name of protein, useful for logging. Uses last
    ::FREESASA_NAME_LIMIT characters if name is too long (since then it
    is probably a file name and the last characters are the most
    interesting). 

    @param s a ::freesasa-object
    @param name Protein name.
*/
    void freesasa_set_proteinname(freesasa *s,const char *name);

/**
    Get protein name. 

    @param s a ::freesasa-object
    @return Protein name.
*/
    const char* freesasa_get_proteinname(const freesasa *s);

/**
    Include or don't include `HETATM` when reading PDB files.

    @param s a ::freesasa-object
    @param include 1 means use `HETATM` records, 0 means don't.
*/
    void freesasa_include_hetatm(freesasa *s, int include);

/////////////
// Results //
/////////////

/**
    Total SASA. 

    Asserts that calculations have been performed.
    
    @param s a ::freesasa-object
    @return Total SASA in Å^2. x
*/
    double freesasa_area_total(const freesasa *s);

/**
    SASA of a certain class of atoms

    Asserts that calculations have been performed and that there is
    information about the atoms (i.e. calculations were not done
    on "abstract" coordinates).

    @param s a ::freesasa-object
    @param c Class of atoms (polar/apolar/nucleic/unknown). 
    @return SASA of class c. 
*/
    double freesasa_area_class(const freesasa *s, freesasa_class c);

/**
    Print SASA for all residue types to file.

    Prints name/value-pairs with the total SASA of each residue
    type. The standard 20 amino acids are always included in output,
    non-standard ones and nucleotides only if they were present in
    input. Each line in the output is prefixed by the string 'RES:'.

    Asserts that calculations have been performed and that there is
    information about the atoms (i.e. calculations were not done
    on "abstract" coordinates).

    @param s a ::freesasa-object
    @param output Output file.
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.
*/
    int freesasa_per_residue_type(const freesasa *s, FILE *output);

/**
    Print SASA for each residue individually to file. 

    Each line in the output is prefixed by the string 'SEQ:'.


    @param s a ::freesasa-object
    @param output Output file.
    @return ::FREESASA_FAIL if problems writing to
    output. ::FREESASA_SUCCESS else.
 */
    int freesasa_per_residue(const freesasa *s, FILE *output);
/**

    Total SASA for specific residue type.

    Returns total SASA for residue of specified type. If the value of
    res_name is unknown, a warning is printed and the SASA value for
    residue type UNK returned, because that is where it would be
    stored. I.e. if residues not known by FreeSASA are used, the only
    option currently is to group them under "UNK". 

    Asserts that calculations have been performed and that there is
    information about the atoms (i.e. calculations were not done
    on "abstract" coordinates).

    @param s a ::freesasa-object
    @param res_name The residue (string of format "ALA", "PHE", etc).
    @return SASA of residue type res_name. 
*/
    double freesasa_area_residue(const freesasa *s, const char *res_name);

/**
    Don't use yet!

    Groups atoms, and returns SASA for each group and a string
    describing it. The available types of atom groups. The returned
    struct should be freed using freesasa_free_strvp().

    @see freesasa_free_strvp()

    @param s a ::freesasa-object
    @param type The type of result
    @return the string-value pairs. Returns NULL if calculation has
    not been performed yet or if type is illegal
 */
    freesasa_strvp* freesasa_string_value_pairs(const freesasa *s,freesasa_result_type type);

/**
    Don't use yet!
    
    Free string-value pair arrays allocated by freesasa_string_value_pairs()

    @param s the object to free.
    @see freesasa_string_value_pairs()
 */
    void freesasa_strvp_free(freesasa_strvp* s);

/**
    Write SASA values to PDB-file.
    
    Takes original PDB and replaces B-factors with those from latest
    calculation. 

    Asserts that calculations have been performed.

    @param s a ::freesasa-object
    @param output File to write to.
    @return ::FREESASA_FAIL if there is no previous PDB input to base
    output on. ::FREESASA_SUCCESS else.
 */
    int freesasa_write_pdb(const freesasa *s, FILE *output);


//////////////////////////////////
// Results for individual atoms //
//////////////////////////////////

/**
    SASA value for given atom.

    Asserts that calculations have been performed and that index is
    valid.

    @param s a ::freesasa-object
    @param i Atom index
    @return SASA value for atom i. 
*/
    double freesasa_area_atom(const freesasa *s, int i);

/**
    SASA for all atoms individually.

    Asserts that calculations have been performed.

    @param s a ::freesasa-object
    @return Array of SASA for all atoms. 
*/
    const double* freesasa_area_atom_array(const freesasa *s);

/**
    Radius of given atom.
    
    @param s a ::freesasa-object
    @param i Atom index.
    @return Radius of atom in Ångström. Prints error and returns negative
    value if atom index is invalid or no value available.
*/
    double freesasa_radius_atom(const freesasa *s, int i);

/**
    Get atomic radii.
    
    @param s a ::freesasa-object
    @return Array of all atomic radii. Returns NULL if no radii are
    available.
*/
    const double* freesasa_radius_atom_array(const freesasa *s);

////////////////////////////
// Other types of results //
////////////////////////////

/**
    Log calculation results.

    Prints log of calculation results to specified file. 

    @param s a ::freesasa-object
    @param log Output-file.
    @return ::FREESASA_SUCCESS on success, ::FREESASA_WARN if
    problems writing to file.
*/
    int freesasa_log(const freesasa *s, FILE *log);

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

#ifdef __cplusplus
}
#endif

#endif
