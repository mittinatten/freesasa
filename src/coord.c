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

struct coord_ {
    double *xyz;
    size_t n;
};

coord_t* coord_new() 
{
    coord_t* c = (coord_t*) malloc(sizeof(coord_t*));
    c->xyz = NULL;
    c->n = 0;
    return c;
}

void coord_free(coord_t *c) 
{
    assert(c != NULL && "NULL-pointer passed to coord_free(1)");
    if (c->xyz) free(c->xyz);
    free(c);
}

coord_t* coord_copy(const coord_t *src) 
{
    assert(src != NULL);
    coord_t *c = coord_new();
    coord_set_all(c,src->xyz,src->n);
    return c;
}

const coord_t* coord_new_linked(double *xyz, size_t n)
{
    assert(xyz != NULL);
    coord_t *c = coord_new();
    c->xyz = xyz;
    c->n = n;
    return c;
}

void coord_append(coord_t *c, const double *xyz, size_t n)
{
    assert(c   != NULL);
    assert(xyz != NULL);

    size_t n_old = c->n;
    c->n += n;
    c->xyz = (double*) realloc(c->xyz, sizeof(double)*3*c->n);
    double *dest = memcpy(&(c->xyz[3*n_old]), xyz, sizeof(double)*n*3);

    assert(dest != NULL);
}

void coord_append_xyz(coord_t *c, const double *x, const double *y, 
		      const double *z, size_t n)
{
    assert(c != NULL);
    assert(x != NULL);
    assert(y != NULL);
    assert(z != NULL);

    double *xyz = (double*) malloc(sizeof(double)*3*n);
    for (int i = 0; i < n; ++i) {
	xyz[i*3] = x[i];
	xyz[i*3+1] = y[i];
	xyz[i*3+2] = z[i];
    }
    coord_append(c,xyz,n);
}

void coord_set_i(coord_t *c, int i, const double* xyz) 
{
    assert(c   != NULL);
    assert(xyz != NULL);
    assert(c->n > i);
    assert(i >= 0);
    
    memcpy(&c->xyz[i*3], xyz, 3*sizeof(double));
}

void coord_set_i_xyz(coord_t *c,int i,double x,double y,double z)
{
    assert(c   != NULL);
    assert(c->n > i);
    assert(i >= 0);

    double *v_i = &c->xyz[i*3];
    *(v_i++) = x;
    *(v_i++) = y;
    *v_i = z;
}

void coord_set_all(coord_t *c,const double* xyz,size_t n) 
{
    assert(c   != NULL);
    assert(xyz != NULL);
    if (c->xyz) free(c->xyz);
    c->n = 0;
    coord_append(c,xyz,n);    
}

void coord_set_all_xyz(coord_t *c,const double* x, const double *y,
		       const double *z, size_t n)
{
    assert(c != NULL);
    if (c->xyz) free(c->xyz);
    c->n = 0;
    coord_append_xyz(c, x, y, z, n);
}

void coord_set_length_i(coord_t *c, int i, double l)
{
    assert(c != NULL);
    assert(c->xyz != NULL);
    
    double x = c->xyz[3*i], y = c->xyz[3*i+1], z = c->xyz[3*i+2];
    double r = sqrt(x*x + y*y + z*z);
    c->xyz[3*i]   *= l/r;
    c->xyz[3*i+1] *= l/r;
    c->xyz[3*i+2] *= l/r;
}

void coord_set_length_all(coord_t *c, double l)
{
    for (int i = 0; i < c->n; ++i) coord_set_length_i(c,i,l);
}

const double* coord_i(const coord_t *c, int i)
{
    assert(c != NULL);
    assert(i < c->n);
    assert(i >= 0);
    return &c->xyz[3*i];
}

double coord_dist(const coord_t *c, int i, int j)
{
    return sqrt(coord_dist2(c,i,j));
}

static inline double dist2(const double *v1, const double *v2)
{
    double dx = v1[0]-v2[0], dy = v1[1]-v2[1], dz = v1[2]-v2[2];
    return dx*dx + dy*dy + dz*dz;    
}

double coord_dist2(const coord_t *c, int i, int j)
{
    double *v1 = &c->xyz[3*i];
    double *v2 = &c->xyz[3*j];
    return dist2(v1,v2);
}

double coord_dist2_12(const coord_t* c1, const coord_t* c2, int i1, int i2)
{
    double *v1 = &c1->xyz[3*i1];
    double *v2 = &c2->xyz[3*i2];
    return dist2(v1,v2);
}

const double* coord_all(const coord_t *c)
{
    assert(c != NULL);
    return c->xyz;
}

size_t coord_n(const coord_t* c)
{
    assert(c != NULL);
    return c->n;
}

void coord_translate(coord_t *c, const double *xyz)
{
    assert(xyz != NULL);
    coord_translate_xyz(c,xyz[0],xyz[1],xyz[2]);
}

void coord_translate_xyz(coord_t *c, double x, double y, double z)
{
    assert(c != NULL);
    
    for (int i = 0; i < c->n; ++i) {
	c->xyz[3*i]   += x; 
	c->xyz[3*i+1] += y; 
	c->xyz[3*i+2] += z; 
    }
}

void coord_scale(coord_t *c, double s)
{
    assert(c != NULL);
    for (int i = 0; i < c->n*3; ++i) {
	c->xyz[i] *= s;
    }
}

