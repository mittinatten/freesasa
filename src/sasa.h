#ifndef SASA_H
#define SASA_H

#include <stdlib.h>
#include "vector3.h"

//don't change this without actually generating new points
#define MAX_SR_POINTS 2000

//this should be made a user option eventually
#define PROBE_RADIUS 1.4 

/** Solvent accessible surface area for each atom is written to the
    array 'sasa'. The user is responsible for making sure this has the
    right size.  'n_points' is the number of points (<= MAX_SR_POINTS)
    to use. Fewer points lead to faster but less accurate
    calculations. */
void sasa_shrake_rupley(double *sasa,
                        const vector3 *xyz,
                        const double *radii,
                        size_t n_atoms,
                        int n_points);

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

#endif
