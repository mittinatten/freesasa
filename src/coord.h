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

/** This is only for internal use, hence all errors are handled with
    asserts */

typedef struct coord_ coord_t;

coord_t* coord_new();

void coord_free(coord_t*);

/** creates a new copy of src */
coord_t* coord_copy(const coord_t *src);

/** Creates a const coord_t-object that is linked to an array of
    coordinates owned by callee. This allows repeated calculations on
    changing object without reinitialization. The const-ness of the
    returned pointer makes sure coordinates can not be modified
    later. It can however not be freed either. */
const coord_t* coord_new_linked(double *xyz, size_t n);

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

void coord_i(double *xyz, int i, const coord_t*);

/** No index bounds checking (for speed) */
double coord_dist(const coord_t*, int i, int j);

/** No index bounds checking (for speed) */
double coord_dist2(const coord_t*, int i, int j);

/** No index bounds checking (for speed) */
int coord_dist2_lt(const coord_t*, int i, int j, double cutoff2); 

const double* coord_all(const coord_t*);

size_t coord_n(const coord_t*);
/*
void coord_add_to(coord_t*, const coord_t*);

coord_t* coord_sum(const coord_t*, const coord_t*);
*/

void coord_translate(coord_t*, const double *xyz);

void coord_translate_xyz(coord_t*, double x, double y, double z);


#endif
