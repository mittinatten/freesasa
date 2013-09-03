/*
  Copyright Simon Mitternacht 2013.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SASA_H
#define SASA_H

#include <stdio.h>
#include "vector3.h"
#include "atomclassifier.h"
#include "protein.h"

//don't change this without actually generating new points
#define MAX_SR_POINTS 2000

//this should be made a user option eventually
#define PROBE_RADIUS 1.4 

/** Solvent accessible surface area for each atom is written to the
    array 'sasa'. The user is responsible for making sure this has the
    right size.  'n_points' is the number of points (<= MAX_SR_POINTS)
    to use. Fewer points lead to faster but less accurate
    calculations. Last argument sets the number of threads in parallel
    computation, only used if the program was compiled with
    -DPTHREADS. */
void sasa_shrake_rupley(double *sasa,
                        const vector3 *xyz,
                        const double *radii,
                        size_t n_atoms,
                        int n_points,
			int n_threads);

/** Solvent accessible surface area for each atom is written to the
    array 'sasa'. The user is responsible for making sure this has the
    right size. The argument grid sets the distance between grid
    points in Å. */
void sasa_lee_richards(double* sasa,
                       const vector3 *xyz,
                       const double *radii,
                       size_t n_atoms,
                       double grid);

// other algorithms?

/** Writes result of SASA_calculation to FILE* summed up for the 
    atom-classes according to the provided atomclassifier. */
void sasa_per_atomclass(FILE*, atomclassifier, 
			protein*, double *sasa);

#endif
