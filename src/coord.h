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

#ifndef FREESASA_COORD_H
#define FREESASA_COORD_H

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

/** This is only for internal use, hence all errors are handled with
    asserts */

typedef struct freesasa_coord_ freesasa_coord_t;

freesasa_coord_t* freesasa_coord_new();

void freesasa_coord_free(freesasa_coord_t*);

/** creates a new copy of src */
freesasa_coord_t* freesasa_coord_copy(const freesasa_coord_t *src);

/** Creates a const freesasa_coord_t-object that is linked to an array of
    coordinates owned by callee. This allows repeated calculations on
    a changing object without reinitialization. The returned
    freesasa_coord_t-pointer is not explicitly const, to allow it to be freed
    later, but objects initiated through this interface will not
    change their coordinates. */
freesasa_coord_t* freesasa_coord_new_linked(const double *xyz, size_t n);

void freesasa_coord_append(freesasa_coord_t*,const double *xyz,size_t n);

void freesasa_coord_append_xyz(freesasa_coord_t*, 
			      const double *x, const double *y, 
			      const double *z, size_t n);

void freesasa_coord_set_i(freesasa_coord_t*,int i,const double* xyz);

void freesasa_coord_set_i_xyz(freesasa_coord_t*,int i,double x,double y,double z);

/* resets everything */
void freesasa_coord_set_all(freesasa_coord_t*,const double* xyz,size_t n);

void freesasa_coord_set_all_xyz(freesasa_coord_t*,
			       const double* x, const double *y,
			       const double *z, size_t n);

void freesasa_coord_set_length_i(freesasa_coord_t*, int i, double l);

void freesasa_coord_set_length_all(freesasa_coord_t *c, double l);

const double* freesasa_coord_i(const freesasa_coord_t*, int i);

/** No index bounds checking (for speed) */
double freesasa_coord_dist(const freesasa_coord_t*, int i, int j) 
    __attrib_pure__;

/** No index bounds checking (for speed) */
double freesasa_coord_dist2(const freesasa_coord_t*, int i, int j) 
    __attrib_pure__;

/** */
double freesasa_coord_dist2_12(const freesasa_coord_t* c1, 
			      const freesasa_coord_t* c2, int i1, int i2)
    __attrib_pure__;

const double* freesasa_coord_all(const freesasa_coord_t*) __attrib_pure__;

size_t freesasa_coord_n(const freesasa_coord_t*) __attrib_pure__;

void freesasa_coord_translate(freesasa_coord_t*, const double *xyz);

void freesasa_coord_translate_xyz(freesasa_coord_t*, 
				 double x, double y, double z);

void freesasa_coord_scale(freesasa_coord_t*, double s);

#undef __attrib_pure__

#endif
