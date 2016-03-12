#ifndef FREESASA_INTERNAL_H
#define FREESASA_INTERNAL_H

#include <stdio.h>
#include "freesasa.h"
#include "coord.h"

/**
    @file
    @author Simon Mitternacht

    Functions to perform the actual SASA calculations.
 */

//! The name of the library, to be used in error messages and logging
extern const char *freesasa_name;

//! Shortcut for memory error generation
#define mem_fail() freesasa_mem_fail(__func__,__FILE__,__LINE__) 

#define fail_msg(msg) freesasa_fail_wloc(__func__,__FILE__,__LINE__,msg)


/**
    Calculate SASA using S&R algorithm.
    
    @param sasa The results are written to this array, the user has to
    make sure it is large enough.
    @param c Coordinates of the object to calculate SASA for.
    @param radii Array of radii for each sphere.
    @param probe_radius Probe radius to be used.
    @param n_points Number of points to be used, must be accepted by 
    freesasa_srp_n_is_valid(). 
    @param n_threads Number of threads to use for parallel computations 
    (only leads to performance improvement for largish objects, see 
    manual). Program has to be compiled with `-DPTHREADS` for this option
    to have any effect.
    
    @return ::FREESASA_SUCCESS on success, ::FREESASA_WARN if multiple
    threads are requested when compiled in single-threaded mode (with
    error message). ::FREESASA_FAIL if memory allocation failure.
*/
int
freesasa_shrake_rupley(double *sasa,
                       const coord_t *c,
                       const double *radii,
                       double probe_radius,
                       int n_points,
                       int n_threads);

/**
    Calculate SASA using L&R algorithm.

    Solvent accessible surface area for each atom is written to the
    array 'sasa'. The user is responsible for making sure this has the
    right size. The argument grid sets the distance between grid
    points in Å. Returns FREESASA_SUCCESS on success, FREESASA_WARN if
    multiple threads are requested when compiled in single-threaded
    mode (with error message). 
    
    @param sasa The results are written to this array, the user has to
    make sure it is large enough.
    @param c Coordinates of the object to calculate SASA for.
    @param radii Array of radii for each sphere.
    @param probe_radius Probe radius to be used.
    @param n_slices_per_atom Number of slices per atom (resolution).
    @param n_threads Number of threads to use for parallel computations 
    (only leads to performance improvement for largish objects, see 
    manual). Program has to be compiled with `-DPTHREADS` for this option
    to have any effect.
    
    @return ::FREESASA_SUCCESS on success, ::FREESASA_WARN if
    multiple threads are requested when compiled in single-threaded
    mode (with error message). ::FREESASA_FAIL if memory allocation 
    failure.
*/
int freesasa_lee_richards(double* sasa,
                          const coord_t *c,
                          const double *radii,
                          double probe_radius,
                          int n_slices_per_atom,
                          int n_threads);


/**
    Get coordinates.
    
    @param s Self.
    @return The coordinates of the structure as a ::coord_t struct.
 */
const coord_t *
freesasa_structure_xyz(const freesasa_structure *s);

/**
    Get number of residues.

    Calculated crudely by determining the number of unique
    combinations of residue name and chain label contained in the
    structure. If residues are mingled i.e. atoms of the same residue
    are in non-contiguous regions of the file, this function might be
    off.

    @param s Self.
    @return Number of residues.
 */
int
freesasa_structure_n_residues(const freesasa_structure *s);

/**
    Get a string describing an atom. 
    Format: "A    1 ALA  CA " 
    (chain label, residue number, residue type, atom name)

    @param s Self.
    @param i Atom index
    @return Descriptor string. 
 */
const char*
freesasa_structure_atom_descriptor(const freesasa_structure *s,
                                   int i);

/**
    Get indices of first and last atoms of a residue
 
    @param s Self.
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
    Get a string describing a residue.
    Format: "A    1 ALA" (chain label, residue number, residue type)
    
    @param s Self.
    @param r_i atom index
    @return Descriptor string
 */
const char*
freesasa_structure_residue_descriptor(const freesasa_structure *s,
                                      int r_i);

/**
    Returns the SASA for a given residue

    @param r The SASA results
    @param s The structure
    @param r_i Index of residue
    @return The SASA of the residue
*/
double
freesasa_single_residue_sasa(const freesasa_result *r,
                             const freesasa_structure *s, 
                             int r_i);
/**
    Holds range in a file, to be initalized with ftell() and used
    with fseek().
 */
struct file_range {
    long begin; //!< Position of beginning of range
    long end; //!< Position of end of range
};

/**
    For convenience, get a file range that covers a whole file.

    @param file The file to study
    @return the ::file_range.
 */
struct file_range
freesasa_whole_file(FILE* file);

/**
    Print failure message using format string and arguments.

    @param format Format string
    @return ::FREESASA_FAIL
*/
int 
freesasa_fail(const char *format,...);

/**
    Print warning message using format string and arguments.

    @param format Format string
    @return ::FREESASA_WARN
 */
int
freesasa_warn(const char *format,...);

/**
    Print warning message using function, file and line-number. 

    Usually used via the macro mem_fail().

    @param func Name of the function that failed
    @param file The file where the function is.
    @param line Line number for the error.
    @return ::FREESASA_FAIL
 */
int
freesasa_mem_fail(const char* func,
                  const char* file,
                  int line);

/**
    Returns string explaining return values of pthread_create() and
    pthread_join().

    @param error_code The error code
    @return A string describing the error code.
 */
const char*
freesasa_thread_error(int error_code);

/**
    Prints fail message with function name, file name, and line number.

    Mainly intended to be used by the macros fail_msg() and mem_fail().

    @param func The function name
    @param file The file name
    @param line The line number
    @param msg The error message
    @return ::FREESASA_FAIL
 */
int
freesasa_fail_wloc(const char* func,
                   const char* file,
                   int line,
                   const char *msg);
                

#endif /* FREESASA_INTERNAL_H */
