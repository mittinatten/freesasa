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

#ifndef SASALIB_COORD_H
#define SASALIB_COORD_H

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

/** This is only for internal use, hence all errors are handled with
    asserts */

typedef struct coord_ coord_t;

coord_t* coord_new();

void coord_free(coord_t*);

/** creates a new copy of src */
coord_t* coord_copy(const coord_t *src);

/** Creates a const coord_t-object that is linked to an array of
    coordinates owned by callee. This allows repeated calculations on
    a changing object without reinitialization. The returned
    coord_t-pointer is not explicitly const, to allow it to be freed
    later, but objects initiated through this interface will not
    change their coordinates. */
coord_t* coord_new_linked(const double *xyz, size_t n);

void coord_append(coord_t*,const double *xyz,size_t n);

void coord_append_xyz(coord_t*, const double *x, const double *y, 
                      const double *z, size_t n);

void coord_set_i(coord_t*,int i,const double* xyz);

void coord_set_i_xyz(coord_t*,int i,double x,double y,double z);

void coord_set_all(coord_t*,const double* xyz,size_t n);

void coord_set_all_xyz(coord_t*,const double* x, const double *y,
                       const double *z, size_t n);

void coord_set_length_i(coord_t*, int i, double l);

void coord_set_length_all(coord_t *c, double l);

const double* coord_i(const coord_t*, int i);

/** No index bounds checking (for speed) */
double coord_dist(const coord_t*, int i, int j) __attrib_pure__;

/** No index bounds checking (for speed) */
double coord_dist2(const coord_t*, int i, int j) __attrib_pure__;

/** */
double coord_dist2_12(const coord_t* c1, const coord_t* c2, int i1, int i2)
    __attrib_pure__;

const double* coord_all(const coord_t*) __attrib_pure__;

size_t coord_n(const coord_t*) __attrib_pure__;

void coord_translate(coord_t*, const double *xyz);

void coord_translate_xyz(coord_t*, double x, double y, double z);

void coord_scale(coord_t*, double s);

#undef __attrib_pure__

#endif
