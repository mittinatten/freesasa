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
 */

//! Struct to store coordinates
typedef struct freesasa_coord freesasa_coord;

//! Initialize new ::freesasa_coord object
freesasa_coord* freesasa_coord_new(void);

//! Free resources allocated by ::freesasa_coord object
void freesasa_coord_free(freesasa_coord*);

/**
   Copy coordinates.

   Creates a new ::freesasa_coord object that is a copy of the
   argument `src`.

   @param src Coordinates to be copied.
   @return Copy of coordinates
*/
freesasa_coord* freesasa_coord_copy(const freesasa_coord *src);

/**
    Creates a `const` ::freesasa_coord-object that is linked to an array of
    coordinates owned by callee. 

    This allows repeated calculations on a changing object without
    reinitialization. The returned ::freesasa_coord-pointer is not
    explicitly const, to allow it to be freed later, but objects
    initiated through this interface will not change their
    coordinates.

    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
    @return New linked ::freesasa_coord object
*/
freesasa_coord* freesasa_coord_new_linked(const double *xyz, size_t n);

/**
    Append coordinates to ::freesasa_coord object from one array.

    @param s Self.
    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
 */
void freesasa_coord_append(freesasa_coord *s,const double *xyz,size_t n);

/**
    Append coordinates to ::freesasa_coord object from three
    separate arrays.

    @param s Self.
    @param x Array of x-coordinates
    @param y Array of x-coordinates
    @param z Array of x-coordinates
    @param n Size of arrays x, y and z.
 */
void freesasa_coord_append_xyz(freesasa_coord *s,
                               const double *x, const double *y,
                               const double *z, size_t n);

/**
    Set given coordinate.

    @param s Self.
    @param i Coordinate to update.
    @param xyz Array with coordinates x,y,z.
 */
void freesasa_coord_set_i(freesasa_coord *s,int i,const double* xyz);

/**
    Set given coordinate.

    @param s Self.
    @param i Coordinate to update.
    @param x x-coordinate.
    @param y y-coordinate.
    @param z z-coordinate
 */
void freesasa_coord_set_i_xyz(freesasa_coord *s,int i,double x,double y,double z);

/**
    Reset everything.
    
    @param s Self.
    @param xyz Array of coordinates x1,y1,z1,x2,y2,z2,...
    @param n Number of coordinates (array has size 3*n).
 */
void freesasa_coord_set_all(freesasa_coord *s,const double* xyz,size_t n);

/**
    Reset everything.
    
    @param s Self.
    @param x Array of x-coordinates
    @param y Array of x-coordinates
    @param z Array of x-coordinates
    @param n Size of arrays x, y and z.
 */
void freesasa_coord_set_all_xyz(freesasa_coord *s,
                                const double* x, const double *y,
                                const double *z, size_t n);

/**
    Set length of a given coordinate vector. Useful for test-points in
    S&R.

    @param s Self.
    @param i The coordinate.
    @param l Desirde length.
 */
void freesasa_coord_set_length_i(freesasa_coord *s, int i, double l);

/**
    Set length of all coordinate vectors to the same length. Useful for test-points in
    S&R.

    @param s Self.
    @param l Desired length.
 */
void freesasa_coord_set_length_all(freesasa_coord *s, double l);

/**
    Coordinates of a given atom.

    @param s Self.
    @param i The coordinate.
    @return Array with coordinates x,y,z.
 */
const double* freesasa_coord_i(const freesasa_coord *s, int i);

/**
    Calculate distance between two coordinate vectors. For speed, the indices i and
    j will not be checked.
    
    @param s Self.
    @param i The first coordinate.
    @param j The second coordinate.
    @return Distance between coorsinate i and j.
*/
double freesasa_coord_dist(const freesasa_coord *s, int i, int j)
    __attrib_pure__;

/**
    Calculate square distance between two coordinate vectors. For
    speed, the indices i and j will not be checked.
    
    @param s Self.
    @param i The first coordinate.
    @param j The second coordinate.
    @return Square distance between coordinate i and j.
*/
double freesasa_coord_dist2(const freesasa_coord *s, int i, int j)
    __attrib_pure__;

/**
    Calculate square distance between two coordinate vectors in
    separate coordinate sets. For speed, the indices i and j will not
    be checked.
    
    @param c1 First set of coordinates..
    @param c2 Second set of coordinates.
    @param i1 First coordinate.
    @param i2 Second coordinate.
    @return Square distance between coordinates i1 and i2.
*/
double freesasa_coord_dist2_12(const freesasa_coord* c1,
                               const freesasa_coord* c2, int i1, int i2)
    __attrib_pure__;

/**
    All coordinates as an array.

    @param s Self.
    @return Array of coordinates x1,y1,z1,x2,y2,z2,...
 */
const double* freesasa_coord_all(const freesasa_coord *s) __attrib_pure__;

/**
    Number of coordinates.
    
    @param s Self.
    @return Number of coordinates.
 */
size_t freesasa_coord_n(const freesasa_coord *s) __attrib_pure__;

/**
    Translate all coordinates by same vector.

    @param s Self.
    @param xyz Array describing translation vector x,y,z.
 */
void freesasa_coord_translate(freesasa_coord *s, const double *xyz);
/**
    Translate all coordinates by same vector.

    @param s Self.
    @param x x-coordinate of translation vector
    @param y y-coordinate of translation vector
    @param z z-coordinate of translation vector
 */
void freesasa_coord_translate_xyz(freesasa_coord *s,
                                  double x, double y, double z);

/**
    Scale all coordinates by given factor.
    
    @param s Self.
    @param a Factor to scale by.
 */
void freesasa_coord_scale(freesasa_coord *s, double a);

#undef __attrib_pure__

#endif
