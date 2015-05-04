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

/** 
    @file
    @author Simon Mitternacht
    
    @section Description

    freesasa.h provides functions to init and perform a SASA
    calculation using FreeSASA with standard atom classifications. The
    user may optionally select algorithm and provide parameters.

    @subsection Coordinates

    If users wish to supply their own coordinates and radii, these are
    accepted as arrays of doubles. The coordinate-array should have
    size 3*n with coordinates in the order x1,y1,z1,x2,y2,z2,...

    @subsection Error-reporting 

    All errors are written to stderr and are prefixed with the string
    'freesasa'. There are two error return values FREESASA_WARN and
    FREESASA_FAIL (see documentation of each function to see when
    these are used). The return value FREESASA_SUCCESS means no errors
    have been spotted.

    Functions that have real-valued return values return negative
    numbers if calculation failed for some reason. The documentation
    for each function explains when this can happen.
*/

#include <stdio.h>

typedef struct freesasa_t freesasa_t;
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
#define FREESASA_FAIL -1 
#define FREESASA_WARN -2 
                         

/*****************************************/
/** Initialization, freeing and copying **/
/*****************************************/

#ifdef __cplusplus
extern "C"{
#endif

/** 
    Initiate freesasa_t-object

    Allocates empty freesasa_t object with default parameters. Must be
    called before calculation. 

    @return The created object
*/
    freesasa_t* freesasa_init();

/** 
    Frees resources allocated to a freesasa_t object. 

    @param s The object to be freed
*/
    void freesasa_free(freesasa_t *s);

/** 
    Copy parameters
    
    Copies all algorithm parameters. Data from calculations are not
    copied, i.e. this can be used for setting up several calculations
    with identical parameters for different proteins.
    
    @param target The target configuration
    @param source The source configuration
*/
    void freesasa_copy_param(freesasa_t *target, const freesasa_t *source);


/**************************/
/** Perform calculations **/
/**************************/

/** 
    Calculate SASA from coordinates and radii. 

    Can be run on any set of spheres. Coordinates and radii will not
    be stored. This will only generate a total SASA, not polar or
    apolar, since only radii are provided, not atom types.

    @see freesasa_calc_pdb()
    @see freesasa_calc_atoms()
    @param s A freesasa_t object. Used for parameters and to store results.
    @param coord An array of coordinates for sphere centers (x1,y1,z1,x2,y2,z2,...)
    @param r An array of radii for the spheres
    @param n number of spheres
    
    @return FREESASA_SUCESS upon successful calculation, prints and
    error and returns FREESASA_FAIL else. 
*/
    int freesasa_calc_coord(freesasa_t  *s, const double *coord,
                            const double *r, size_t n);

/** 
    Calculate SASA from PDB-file.
    
    Reads PDB-file and calculates SASA. HETATM records are
    ignored. Results stored in parameter s. If s is not initialized
    default values are used, these are stored in s.  If the object has
    been used in calculations previously, the results from these will
    be over-written. 

    @see freesasa_calc_atoms()
    @see freesasa_calc_coord()
    @param s A freesasa_t object. Used for parameters and to store results.
    @param pdb_file PDB-file for input
    @return FREESASA_SUCCESS if calculation successful, prints an
    error and returns FREESASA_FAIL if not.
*/
    int freesasa_calc_pdb(freesasa_t *s, FILE *pdb_file);

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
    @param s A freesasa_t object. Used for parameters and to store results.
    @param coord Array of atom coordinates (x1,y1,z1,x2,y2,z2,...). Should
      have 3*n elements.
    @param resnames Array of strings of the format "ALA", "PHE", etc. 
      Should have n elements.
    @param atomnames Array of strings of the format " CA ", " OXT", etc. 
      Should have n elements.
    @param number of atoms
    @return FREESASA_SUCCESS if calculation successful.  FREESASA_WARN
      if atoms or coordinates have invalid formats or if any atoms can't
      be classified. FREESASA_FAIL if calculations failed.
*/
    int freesasa_calc_atoms(freesasa_t *s, const double *coord,
                            const char **resnames, 
                            const char **atomnames, size_t n);

/** Reads pdb-file and calculates radii for each atom. Memory is
    allocated to store them in the array 'r'. The return value is the
    size of the allocated array. Prints error and returns FREESASA_FAIL
    if reading input fails. This can be used if coordinates are linked
    as below and default radii are to be used.
    Not properly tested yet!
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
    Link a set of coordinates to the freesasa_t object.

    If the linked coordinates are updated freesasa_refresh() can be
    used to recalculate SASA. FreeSASA will not change the
    coordinates. If the freesasa_t-object has been initalized with a
    PDB file the corresponding memory is released and the results
    erased.
    
    @see freesasa_refresh()
    @param s A freesasa_t object. Used for parameters and to store results.
    @param coord An array of coordinates for sphere centers (x1,y1,z1,x2,y2,z2,...)
    @param r An array of radii for the spheres
    @param n number of spheres
    @return At the moment always returns FREESASA_SUCCESS
*/
    int freesasa_link_coord(freesasa_t *s, const double *coord,
                            double *r, size_t n);

/** 
    Recalculate SASA based on linked coordinates.

    Recalculates SASA, based on the assumption that a set of external
    coordinates have been updated elsewhere. 

    @see freesasa_link_coord()
    @param s A freesasa_t object. Used for parameters, coordinate-link
    and to store results.  
    @return FREESASA_FAIL if no coordinates or radii are found in
    s. FREESASA_SUCCESS upon successful computation.
*/
    int freesasa_refresh(freesasa_t *s);

/** 
    The number of atoms in the latest SASA calculation.  
    
    Returns 0 if no coordinates have been linked or no calculation
    has been performed.

    @param s Self.
    @return The number of atoms.
*/
    size_t freesasa_n_atoms(const freesasa_t *s);

/******************************/
/** Settings for calculation **/
/******************************/

/** 
    Set algorithm. 

    @param s Self.
    @param alg The algorithm to be set.
    @return FREESASA_SUCCESS if alg is valid, FREESASA_WARN else. 
*/
    int freesasa_set_algorithm(freesasa_t *s, freesasa_algorithm alg);

/** 
    Get algorithm. 

    @param s Self.
    @return The algorithm.
*/
    freesasa_algorithm freesasa_get_algorithm(const freesasa_t *s);

/** 
    Get algorithm as string.
 
    @param s Self.
    @return The algorithm as a string.
*/
    const char* freesasa_algorithm_name(const freesasa_t *s);

/** 
    Set probe radius.

    Sets probe radius for SASA calculations (default FREESASA_DEF_PROBE_RADIUS = 1.4 Å). If submitted radius
    is invalid, default is used and an error message printed.

    @param s Self.
    @param r Value for probe radius in Ångström.
    @return FREESASA_SUCCESS for valid r-values. FREESASA_WARN else. 
*/
    int freesasa_set_probe_radius(freesasa_t *s,double r);

/** 
    Get probe radius.

    @param s Self.
    @return Probe radius in Ångström. 
*/
    double freesasa_get_probe_radius(const freesasa_t *s);
    
/**
    Set number of test points for S&R algorithm.

    Sets number of test points for S&R algorithm (default
    FREESASA_DEF_SR_N = 100).
    
    @param s Self.
    @param n Number of points.
    @return FREESASA_SUCCESS if n is valid. Else: prints error
    message, returns FREESASA_WARN and sets to default value.
*/
    int freesasa_set_sr_points(freesasa_t *s, int n);

/** 
    Get number of test points for S&R algorithm. 

    @param s Self.
    @return Number of points. FREESASA_WARN if S&R algorithm not
    selected.
*/
    int freesasa_get_sr_points(const freesasa_t *s);

/** 
    Set L&R slice width.

    Sets slice width d for L&R algorithm in Ångström (default
    FREESASA_DEF_LR_D = 0.25 Å). 

    @param s Self.
    @param d Slice width.
    @return FREESASA_SUCCESS if d is valid. Else: prints error
    message, returns FREESASA_WARN and sets to default value
    
*/
    int freesasa_set_lr_delta(freesasa_t *s, double d);

/** 
    Get L&R slice width.
    
    Get slice width for L&R algorithm in Ångström. 

    @param s Self.
    @return Slice width. Negative value if L&R algorithm not
    selected. 
*/
    double freesasa_get_lr_delta(const freesasa_t *s);

/** 
    Set number of threads.

    Sets the number of threads to use for parallel computations,
    useful for large proteins and high resolution. Minimal gain for
    smaller systems.

    @param s Self.
    @param n Number of threads.
    @return FREESASA_SUCCESS if n is valid, uses default value and
    returns FREESASA_WARN else.
*/
    int freesasa_set_nthreads(freesasa_t *s,int n);

/** 
    Get number of threads.

    Gets the number of threads used in the calcultion.

    @param s Self.
    @return Number of threads.
*/
    int freesasa_get_nthreads(const freesasa_t *s);

/** 
    Set protein name
    
    Sets name of protein, useful for logging. Uses last
    FREESASA_NAME_LIMIT characters if name is too long (since then it
    is probably a file name and the last characters are the most
    interesting). 

    @param s Self.
    @param name Protein name.
*/
    void freesasa_set_proteinname(freesasa_t *s,const char *name);

/** 
    Get protein name. 

    @param s Self.
    @return Protein name.
*/
    const char* freesasa_get_proteinname(const freesasa_t *s);

/*************/
/** Results **/
/*************/

/** 
    Total SASA
    
    @param s Self.
    @return Total SASA in Å^2. Negative return value and warning printed
    if calculation hasn't been performed yet. 
*/
    double freesasa_area_total(const freesasa_t *s);

/** 
    SASA of a certain class of atoms

    @param s The freesasa_t object
    @param c Class of atoms (polar/apolar/nucleic/unknown). 
    @return SASA of class c. Negative value if calculation has not
    been performed yet. 
*/
    double freesasa_area_class(const freesasa_t *s, freesasa_class c);

/** 
    Print SASA for all residue types to file.

    Prints name/value-pairs with the total SASA of each residue
    type. The standard 20 amino acids are always included in output,
    non-standard ones and nucleotides only if they were present in
    input. 
    
    @param output Output file.
    @param s Self.

    @return FREESASA_FAIL if file-pointer NULL or if no calculation
    has been performed yet.  
*/
    int freesasa_per_residue(FILE *output, const freesasa_t *s);

/** 
    Total SASA for specific residue type.

    Returns total SASA for residue of specified type. If the value of
    res_name is unknown, a warning is printed and the SASA value for
    residue type UNK returned, because that is where it would be
    stored. I.e. if residues not known by FreeSASA are used, the only
    option currently is to group them under "UNK". 

    @param s Self.
    @param res_name The residue (string of format "ALA", "PHE", etc).

    @return SASA of residue type res_name. Negative value if called
    before calculations have been performed.
*/
    double freesasa_area_residue(const freesasa_t *s, const char *res_name);

/** 
    Write SASA values to PDB-file.
    
    Takes original PDB and replaces B-factors with those from latest
    calculation. 

    @param output File to write to.
    @param s Self.
    @return FREESASA_FAIL if there is no previous PDB input to base
    output on, if there are problems with the output destination, if
    there are no SASA-values to use, or there are inconsistencies
    between stored structure and SASA-values. FREESASA_SUCCESS else.
 */
    int freesasa_write_pdb(FILE *output, const freesasa_t *s);

/**********************************/
/** Results for individual atoms **/
/**********************************/

/** 
    SASA value for given atom.

    @param s Self.
    @param i Atom index
    @return SASA value for atom i. Prints error and returns negative
    value if atom index is invalid or if no calculation has been
    performed. 
*/
    double freesasa_area_atom(const freesasa_t *s, int i);

/** 
    SASA for all atoms individually.
    
    @param s Self.
    @return Array of SASA for all atoms. Returns NULL if no results
    available.
*/
    const double* freesasa_area_atom_array(const freesasa_t *s);

/** 
    Radius of given atom.
    
    @param s Self.
    @param i Atom index.
    @return Radius of atom in Ångström. Prints error and returns negative
    value if atom index is invalid or no value available. 
*/
    double freesasa_radius_atom(const freesasa_t *s, int i);

/** 
    Get atomic radii.
    
    @param s Self.
    @return Array of all atomic radii. Returns NULL if no radii are
    available.
*/
    const double* freesasa_radius_atom_array(const freesasa_t *s);

/****************************/
/** Other types of results **/
/****************************/

/** 
    Log calculation results.

    Prints log of calculation results to specified file. 

    @param s Self.
    @return FREESASA_SUCCESS on success, FREESASA_WARN if
    inconsistencies are detected (with explanatory error-message). 
*/
    int freesasa_log(FILE *log, const freesasa_t *s);

#ifdef __cplusplus
}
#endif

#endif
