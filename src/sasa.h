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

#ifndef FREESASA_SASA_H
#define FREESASA_SASA_H

#include <stdio.h>
#include "coord.h"

/**
    @file
    @author Simon Mitternacht

    Functions to perform the actual SASA calculations.
 */

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


#endif
