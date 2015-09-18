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

#ifndef FREESASA_COORD_H
#define FREESASA_COORD_H

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

/**
    @file
    @author Simon Mitternacht

    This is only for interal use, error handling is done by asserts.
    The header provides functions to store, copy and modify
    coordinates through the type ::freesasa_coord.
    
    The distance calculation functions (freesasa_dist(),
    freesasa_dist2() and freesasa_dist2_12()) are useful for code that
    is not performance critical. When efficieny is a priority it is
    better to use freesasa_coord_all() to obtain a pointer to the
    array of coordinates and calculate directly using this..
 */

//! Struct to store coordinates
typedef struct freesasa_coord freesasa_coord;

/**
    Initialize new ::freesasa_coord object.

    Return value is dynamically allocated, should be freed with
    freesasa_coord_free().
    
    @return An empty ::freesasa_coord object. Returns NULL if out of
    memory.
 */
freesasa_coord * 
freesasa_coord_new(void);

/**
   Free resources allocated by ::freesasa_coord object.

   Will not free the coordinate array itself if it was initialized by
   freesasa_coord_new_linked().
 */
void
freesasa_coord_free(freesasa_coord *coord);

/**
   Copy coordinates.

   Creates a new ::freesasa_coord object that is a copy of the
   argument `src`.

   Return value is dynamically allocated, should be freed with
   freesasa_coord_free().

   @param src Coordinates to be copied.
   @return Copy of coordinates. NULL if out of memory.
 */
freesasa_coord *
freesasa_coord_copy(const freesasa_coord *src);

/**
    Creates a `const` ::freesasa_coord-object that is linked to an array of
    coordinates owned by callee.

    This allows repeated calculations on a changing object without
    reinitialization. The returned ::freesasa_coord-pointer is not
    explicitly const, to allow it to be freed later, but objects
    initiated through this interface will not change their
    coordinates.

    Return value is dynamically allocated, should be freed with
    freesasa_coord_free().

    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return New linked ::freesasa_coord object. NULL if out of memory.
 */
freesasa_coord *
freesasa_coord_new_linked(const double *xyz,
                          int n);

/**
    Append coordinates to ::freesasa_coord object from one array.

    @param coord A ::freesasa_coord object
    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int
freesasa_coord_append(freesasa_coord *coord,
                      const double *xyz,
                      int n);

/**
    Append coordinates to ::freesasa_coord object from three
    separate arrays.

    @param coord A ::freesasa_coord object
    @param x Array of x-coordinates
    @param y Array of x-coordinates
    @param z Array of x-coordinates
    @param n Size of arrays x, y and z.
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int
freesasa_coord_append_xyz(freesasa_coord *coord,
                          const double *x,
                          const double *y,
                          const double *z,
                          int n);

/**
    Set given coordinate.

    Bounds-checking of i is handled by asserts for efficiency,
    i.e. only done in debug-mode.

    @param coord A ::freesasa_coord object
    @param i Index
    @param xyz Array with coordinates x,y,z.
 */
void
freesasa_coord_set_i(freesasa_coord *coord,
                     int i,
                     const double* xyz);

/**
    Set given coordinate.

    Bounds-checking of i is handled by asserts for efficiency,
    i.e. only done in debug-mode.

    @param coord A ::freesasa_coord object
    @param i Index
    @param x x-coordinate.
    @param y y-coordinate.
    @param z z-coordinate
 */
void
freesasa_coord_set_i_xyz(freesasa_coord *coord,
                         int i,
                         double x,
                         double y,
                         double z);

/**
    Reset everything.

    Allocates memory to allow array of size n.
    
    @param coord A ::freesasa_coord object
    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int
freesasa_coord_set_all(freesasa_coord *coord,
                       const double* xyz,
                       int n);

/**
    Reset everything.

    Allocates memory to allow array of size n.    

    @param coord A ::freesasa_coord object
    @param x Array of x-coordinates
    @param y Array of x-coordinates
    @param z Array of x-coordinates
    @param n Size of arrays x, y and z.
    @return FREESASA_SUCCESS if successful, FREESASA_FAIL if out of memory.
 */
int
freesasa_coord_set_all_xyz(freesasa_coord *coord,
                           const double* x,
                           const double *y,
                           const double *z,
                           int n);

/**
    Set length of a given coordinate vector. Useful for test-points in
    S&R.

    Bounds-checking of `i` and `length` is handled by asserts for
    efficiency, i.e. only done in debug-mode.

    @param coord A ::freesasa_coord object
    @param i Index
    @param length Desired length (>= 0)
 */
void
freesasa_coord_set_length_i(freesasa_coord *coord,
                            int i,
                            double length);

/**
    Set length of all coordinate vectors to the same length. 

    This means all coordinates are on the same sphere.

    Bounds-checking of length is handled by asserts for
    efficiency, i.e. only done in debug-mode.

    @param coord A ::freesasa_coord object
    @param length Desired length (>= 0)
 */
void
freesasa_coord_set_length_all(freesasa_coord *coord,
                              double length);

/**
    Coordinates for a given index.

    Bounds-checking of `i` is handled by asserts for efficiency,
    i.e. only done in debug-mode.

    @param coord A ::freesasa_coord object
    @param i Index
    @return Array with coordinates x,y,z
 */
const double*
freesasa_coord_i(const freesasa_coord *coord,
                 int i);

/**
    Calculate distance between two coordinate vectors. For speed,
    arguments aren't checked.
    
    @param coord A ::freesasa_coord object
    @param i first index
    @param j second index
    @return Distance between coorsinate i and j.
*/
double freesasa_coord_dist(const freesasa_coord *coord,
                           int i,
                           int j)
    __attrib_pure__;

/**
    Calculate square distance between two coordinate vectors. For
    speed, arguments aren't checked.
    
    @param coord A ::freesasa_coord object
    @param i First index
    @param j Second index
    @return Square distance between coordinate i and j
*/
double
freesasa_coord_dist2(const freesasa_coord *coord,
                     int i,
                     int j)
    __attrib_pure__;

/**
    Calculate square distance between two coordinate vectors in
    separate coordinate sets. For speed, arguments aren't checked.
    
    @param c1 First set of coordinates
    @param c2 Second set of coordinates
    @param i1 Index in first set
    @param i2 Index in second set
    @return Square distance between coordinates i1 and i2
*/
double freesasa_coord_dist2_12(const freesasa_coord* c1,
                               const freesasa_coord* c2,
                               int i1,
                               int i2)
    __attrib_pure__;

/**
    All coordinates as an array.

    @param coord A ::freesasa_coord object
    @return Array of coordinates x1,y1,z1,x2,y2,z2,...
 */
const double*
freesasa_coord_all(const freesasa_coord *coord) __attrib_pure__;

/**
    Number of coordinates.
    
    @param coord A ::freesasa_coord object
    @return Number of coordinates
 */
int
freesasa_coord_n(const freesasa_coord *coord) __attrib_pure__;

/**
    Translate all coordinates by same vector.

    @param coord A ::freesasa_coord object
    @param xyz Array describing translation vector x,y,z.
 */
void
freesasa_coord_translate(freesasa_coord *coord,
                         const double *xyz);

/**
    Translate all coordinates by same vector.

    @param coord A ::freesasa_coord object
    @param x x-coordinate of translation vector
    @param y y-coordinate of translation vector
    @param z z-coordinate of translation vector
 */
void
freesasa_coord_translate_xyz(freesasa_coord *coord,
                             double x,
                             double y,
                             double z);

/**
    Scale all coordinates by given factor.
    
    @param coord A ::freesasa_coord object
    @param a Factor to scale by
 */
void
freesasa_coord_scale(freesasa_coord *coord,
                     double a);

#undef __attrib_pure__

#endif
