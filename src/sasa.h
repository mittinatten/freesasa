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

#ifndef FREESASA_SASA_H
#define FREESASA_SASA_H

#include <stdio.h>
#include "coord.h"

/** Solvent accessible surface area for each atom is written to the
    array 'sasa'. The user is responsible for making sure this has the
    right size.  'n_points' is the number of points (<= MAX_SR_POINTS)
    to use. Fewer points lead to faster but less accurate
    calculations. Last argument sets the number of threads in parallel
    computation, only used if the program was compiled with
    -DPTHREADS. Returns FREESASA_SUCCESS on success, FREESASA_WARN if
    multiple threads are requested when compiled in single-threaded
    mode (with error message). */
int freesasa_shrake_rupley(double *sasa,
                           const freesasa_coord_t *c,
                           const double *radii,
                           double probe_radius,
                           int n_points,
                           int n_threads);

/** Solvent accessible surface area for each atom is written to the
    array 'sasa'. The user is responsible for making sure this has the
    right size. The argument grid sets the distance between grid
    points in Å. Returns FREESASA_SUCCESS on success, FREESASA_WARN if
    multiple threads are requested when compiled in single-threaded
    mode (with error message). */
int freesasa_lee_richards(double* sasa,
                          const freesasa_coord_t *c,
                          const double *radii,
                          double probe_radius,
                          double grid,
                          int n_threads);

// other algorithms?

#endif
