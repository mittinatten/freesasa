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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "coord.h"

struct sasalib_coord_ {
    double *xyz;
    size_t n;
    int is_const;
};

sasalib_coord_t* sasalib_coord_new() 
{
    sasalib_coord_t* c = (sasalib_coord_t*) malloc(sizeof(sasalib_coord_t));
    c->xyz = NULL;
    c->n = 0;
    c->is_const = 0;
    return c;
}

void sasalib_coord_free(sasalib_coord_t *c) 
{
    assert(c != NULL && "NULL-pointer passed to sasalib_coord_free(1)");
    if (c->xyz && !c->is_const) free(c->xyz);
    free(c);
}

sasalib_coord_t* sasalib_coord_copy(const sasalib_coord_t *src) 
{
    assert(src != NULL);
    sasalib_coord_t *c = sasalib_coord_new();
    sasalib_coord_set_all(c,src->xyz,src->n);
    return c;
}

sasalib_coord_t* sasalib_coord_new_linked(const double *xyz, size_t n)
{
    assert(xyz != NULL);
    sasalib_coord_t *c = sasalib_coord_new();
    c->xyz = (double*)xyz; 
    c->n = n;
    c->is_const = 1;
    return c;
}

void sasalib_coord_append(sasalib_coord_t *c, const double *xyz, size_t n)
{
    assert(c   != NULL);
    assert(xyz != NULL);
    assert(!c->is_const);

    size_t n_old = c->n;
    c->n += n;
    c->xyz = (double*) realloc(c->xyz, sizeof(double)*3*c->n);

    double *dest = memcpy(&(c->xyz[3*n_old]), xyz, sizeof(double)*n*3);

    assert(dest != NULL);
}

void sasalib_coord_append_xyz(sasalib_coord_t *c, 
			      const double *x, const double *y, 
			      const double *z, size_t n)
{
    assert(c != NULL);
    assert(x != NULL);
    assert(y != NULL);
    assert(z != NULL);
    assert(!c->is_const);
    
    double *xyz = (double*)malloc(sizeof(double)*n*3);
    for (int i = 0; i < n; ++i) {
        xyz[i*3] = x[i];
        xyz[i*3+1] = y[i];
        xyz[i*3+2] = z[i];
    }
    sasalib_coord_append(c,xyz,n);
    free(xyz);
}

void sasalib_coord_set_i(sasalib_coord_t *c, int i, const double* xyz) 
{
    assert(c   != NULL);
    assert(xyz != NULL);
    assert(c->n > i);
    assert(i >= 0);
    assert(!c->is_const);
    
    memcpy(&c->xyz[i*3], xyz, 3*sizeof(double));
}

void sasalib_coord_set_i_xyz(sasalib_coord_t *c,int i,
			     double x,double y,double z)
{
    assert(c   != NULL);
    assert(c->n > i);
    assert(i >= 0);
    assert(!c->is_const);

    double *v_i = &c->xyz[i*3];
    *(v_i++) = x;
    *(v_i++) = y;
    *v_i = z;
}

void sasalib_coord_set_all(sasalib_coord_t *c, const double* xyz, size_t n) 
{
    assert(c   != NULL);
    assert(xyz != NULL);
    assert(!c->is_const);

    if (c->xyz) free(c->xyz);
    c->n = 0;
    sasalib_coord_append(c,xyz,n);    
}

void sasalib_coord_set_all_xyz(sasalib_coord_t *c,
			       const double* x, const double *y,
			       const double *z, size_t n)
{
    assert(c != NULL);
    assert(!c->is_const);
    if (c->xyz) free(c->xyz);
    c->n = 0;
    sasalib_coord_append_xyz(c, x, y, z, n);
}

void sasalib_coord_set_length_i(sasalib_coord_t *c, int i, double l)
{
    assert(c != NULL);
    assert(c->xyz != NULL);
    assert(!c->is_const);
    
    double x = c->xyz[3*i], y = c->xyz[3*i+1], z = c->xyz[3*i+2];
    double r = sqrt(x*x + y*y + z*z);
    c->xyz[3*i]   *= l/r;
    c->xyz[3*i+1] *= l/r;
    c->xyz[3*i+2] *= l/r;
}

void sasalib_coord_set_length_all(sasalib_coord_t *c, double l)
{
    assert(c != NULL);
    assert(!c->is_const);
    for (int i = 0; i < c->n; ++i) sasalib_coord_set_length_i(c,i,l);
}

const double* sasalib_coord_i(const sasalib_coord_t *c, int i)
{
    assert(c != NULL);
    assert(i < c->n);
    assert(i >= 0);
    return &c->xyz[3*i];
}

double sasalib_coord_dist(const sasalib_coord_t *c, int i, int j)
{
    return sqrt(sasalib_coord_dist2(c,i,j));
}

static inline double dist2(const double *v1, const double *v2)
{
    double dx = v1[0]-v2[0], dy = v1[1]-v2[1], dz = v1[2]-v2[2];
    return dx*dx + dy*dy + dz*dz;    
}

double sasalib_coord_dist2(const sasalib_coord_t *c, int i, int j)
{
    double *v1 = &c->xyz[3*i];
    double *v2 = &c->xyz[3*j];
    return dist2(v1,v2);
}

double sasalib_coord_dist2_12(const sasalib_coord_t* c1, 
			      const sasalib_coord_t* c2, int i1, int i2)
{
    double *v1 = &c1->xyz[3*i1];
    double *v2 = &c2->xyz[3*i2];
    return dist2(v1,v2);
}

const double* sasalib_coord_all(const sasalib_coord_t *c)
{
    assert(c != NULL);
    return c->xyz;
}

size_t sasalib_coord_n(const sasalib_coord_t* c)
{
    assert(c != NULL);
    return c->n;
}

void sasalib_coord_translate(sasalib_coord_t *c, const double *xyz)
{
    assert(!c->is_const);
    assert(xyz != NULL);
    sasalib_coord_translate_xyz(c,xyz[0],xyz[1],xyz[2]);
}

void sasalib_coord_translate_xyz(sasalib_coord_t *c, 
				 double x, double y, double z)
{
    assert(c != NULL);
    assert(!c->is_const);
    
    for (int i = 0; i < c->n; ++i) {
        c->xyz[3*i]   += x; 
        c->xyz[3*i+1] += y; 
        c->xyz[3*i+2] += z; 
    }
}

void sasalib_coord_scale(sasalib_coord_t *c, double s)
{
    assert(c != NULL);
    assert(!c->is_const);
    for (int i = 0; i < c->n*3; ++i) {
        c->xyz[i] *= s;
    }
}

