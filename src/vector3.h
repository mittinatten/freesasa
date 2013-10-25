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

#ifndef SASALIB_VECTOR3_H
#define SASALIB_VECTOR3_H

#include <math.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

typedef struct {
    double x;
    double y;
    double z;
} vector3;

extern const vector3 vector_zero;
extern const vector3 vector_ex;
extern const vector3 vector_ey;
extern const vector3 vector_ez;

/** set coordinates of a 3-vector */
void vector3_set(vector3* v, double x, double y, double z);

/** set cartesian coordinates based on polar coordinates */
void vector3_set_polar(vector3* v, double polar, double az, double r);

/** length of 3-vector */
double vector3_mag(const vector3* v);

/** length squared of 3-vector */
double vector3_mag2(const vector3* v);

/** multiply length of vector */
void vector3_multiply(vector3* v, double m);

/** set length of vector */
void vector3_setlength(vector3* v, double L);

/** calculate distance-vector v2-v1 */
void vector3_diff(vector3* diff, const vector3* v1, const vector3* v2);

/** add two 3-vectors */
void vector3_sum(vector3* sum, const vector3* v1, const vector3* v2);

/** add v2 to v1 */
void vector3_add(vector3* v1, const vector3* v2);

/** subtract v2 from v1 */
void vector3_sub(vector3* v1, const vector3* v2);

/** calculates the distance squared between two positions */
double vector3_dist2(const vector3* v1, const vector3* v2);

/** calculates the distance squared between two positions over
    periodic boundary conditions (box with side L) */
double vector3_perdist2(const vector3*, const vector3*, const double L);

/** create a copy of a 3-vector*/
void vector3_copy(vector3* u, const vector3* v);

/** scalar product */
double vector3_dot(const vector3* v1, const vector3* v2);

/** vector-product */
void vector3_cross(vector3* u, const vector3* v1, const vector3* v2);

/** angle between v1 and v2 */
double vector3_theta(const vector3* v1, const vector3* v2);

/** measure rotation of v3 around v2 compared to v1 */
double vector3_torsion(const vector3* v1, const vector3* v2, const vector3* v3);

/** rotate vectors in v angle phi around axis at origin from index i1
    to i2 (inclusive) */
void vector3_rotate_block(vector3* v, int i1, int i2,
                          const vector3* axis,
                          const vector3* origin,
                          double phi);

/** For systems with symmetric periodic boundary conditions, how much
    must be added to v to put it in box of side L. Adding the
    calculated vector diff to v will achieve this. */
void vector3_distance_to_box(vector3* diff, const vector3* v, double L);

/** Returns the closest distance between p3 and a line passing p1 and
    p3*/
double vector3_passage_distance(const vector3 *p1,
                                const vector3 *p2,
                                const vector3 *p3);


/** Allocates an array of vector3 from arrays of x, y and z
    coordinates. Array can be freed with regular free(). */
vector3* vector3_xyz_alloc(const double *x, const double *y, 
			   const double *z, int n);

/** Allocates an array of vector3 from an array of size 3*n with
    coordinates in order x1,y1,z1,x2,y2,z2. Array can be freed with
    regular free(). */
vector3* vector3_x3n_alloc(const double *x3n, int n);

#endif
