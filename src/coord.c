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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "freesasa.h"
#include "coord.h"
#include "util.h"

struct freesasa_coord {
    double *xyz;
    int n;
    int is_const;
};

freesasa_coord* freesasa_coord_new()
{
    freesasa_coord* c = malloc(sizeof(freesasa_coord));
    if (c != NULL) {
        c->xyz = NULL;
        c->n = 0;
        c->is_const = 0;
    } else {
        mem_fail();
    }
    return c;
}

void freesasa_coord_free(freesasa_coord *c)
{
    assert(c);
    if (c->xyz && !c->is_const) free(c->xyz);
    free(c);
}
static void coord_clear(freesasa_coord *c)
{
    assert(c);
    assert(!c->is_const);
    if (c->xyz) {
        free(c->xyz);
        c->xyz = NULL;
    }
    c->n = 0;
}


freesasa_coord* freesasa_coord_copy(const freesasa_coord *src)
{
    assert(src);
    freesasa_coord *c = freesasa_coord_new();
    if (c != NULL) freesasa_coord_set_all(c,src->xyz,src->n);
    else mem_fail();
    return c;
}

freesasa_coord* freesasa_coord_new_linked(const double *xyz, int n)
{
    assert(xyz);
    assert(n > 0);
    freesasa_coord *c = freesasa_coord_new();
    if (c != NULL) {
        c->xyz = (double*)xyz;
        c->n = n;
        c->is_const = 1;
    } else mem_fail();
    return c;
}

int freesasa_coord_append(freesasa_coord *c, const double *xyz, int n)
{
    assert(c); assert(xyz); assert(!c->is_const);

    double *dest;
    int n_old = c->n;

    if (n == 0) return FREESASA_SUCCESS;

    c->n += n;
    c->xyz = (double*) realloc(c->xyz, sizeof(double)*3*c->n);
    if (c->xyz == NULL) return mem_fail();

    dest = memcpy(&(c->xyz[3*n_old]), xyz, sizeof(double)*n*3);
    return FREESASA_SUCCESS;
}

int freesasa_coord_append_xyz(freesasa_coord *c,
                               const double *x, const double *y,
                               const double *z, int n)
{
    assert(c); assert(x); assert(y); assert(z);
    assert(!c->is_const);
    double *xyz;
    
    if (n == 0) return FREESASA_SUCCESS;
    
    xyz = malloc(sizeof(double)*n*3);
    if(xyz == NULL) return mem_fail();
    
    for (int i = 0; i < n; ++i) {
        xyz[i*3] = x[i];
        xyz[i*3+1] = y[i];
        xyz[i*3+2] = z[i];
    }
    if (freesasa_coord_append(c,xyz,n) == FREESASA_SUCCESS) {
        free(xyz);
        return FREESASA_SUCCESS;
    }
    return mem_fail();
}

void freesasa_coord_set_i(freesasa_coord *c, int i, const double* xyz)
{
    assert(c); assert(xyz);
    assert(i < c->n && i >= 0);
    assert(!c->is_const);

    memcpy(&c->xyz[i*3], xyz, 3*sizeof(double));
}

void freesasa_coord_set_i_xyz(freesasa_coord *c,int i,
                              double x,double y,double z)
{
    assert(c); assert(c->n > i); assert(i >= 0);
    assert(!c->is_const);

    double *v_i = &c->xyz[i*3];
    *(v_i++) = x;
    *(v_i++) = y;
    *v_i = z;
}

void freesasa_coord_set_all(freesasa_coord *c, const double* xyz, int n)
{
    assert(c); assert(xyz);
    coord_clear(c);
    freesasa_coord_append(c,xyz,n);
}

void freesasa_coord_set_all_xyz(freesasa_coord *c,
                                const double* x, const double *y,
                                const double *z, int n)
{
    assert(c); assert(x); assert(y); assert(z);
    coord_clear(c);
    freesasa_coord_append_xyz(c, x, y, z, n);
}

void freesasa_coord_set_length_i(freesasa_coord *c, int i, double l)
{
    assert(c); assert(c->xyz);
    assert(!c->is_const);
    assert(i >= 0 && i < c->n);
    assert(l >= 0);

    double x = c->xyz[3*i], y = c->xyz[3*i+1], z = c->xyz[3*i+2];
    double r = sqrt(x*x + y*y + z*z);
    c->xyz[3*i]   *= l/r;
    c->xyz[3*i+1] *= l/r;
    c->xyz[3*i+2] *= l/r;
}

void freesasa_coord_set_length_all(freesasa_coord *c, double l)
{
    assert(c);
    assert(!c->is_const);
    for (int i = 0; i < c->n; ++i) freesasa_coord_set_length_i(c,i,l);
}

const double* freesasa_coord_i(const freesasa_coord *c, int i)
{
    assert(c);
    assert(i < c->n);
    assert(i >= 0);
    return &c->xyz[3*i];
}

double freesasa_coord_dist(const freesasa_coord *c, int i, int j)
{
    return sqrt(freesasa_coord_dist2(c,i,j));
}

static inline double dist2(const double *v1, const double *v2)
{
    double dx = v1[0]-v2[0], dy = v1[1]-v2[1], dz = v1[2]-v2[2];
    return dx*dx + dy*dy + dz*dz;
}

double freesasa_coord_dist2(const freesasa_coord *c, int i, int j)
{
    double *v1 = &c->xyz[3*i];
    double *v2 = &c->xyz[3*j];
    return dist2(v1,v2);
}

double freesasa_coord_dist2_12(const freesasa_coord* c1,
                               const freesasa_coord* c2, int i1, int i2)
{
    double *v1 = &c1->xyz[3*i1];
    double *v2 = &c2->xyz[3*i2];
    return dist2(v1,v2);
}

const double* freesasa_coord_all(const freesasa_coord *c)
{
    assert(c);
    return c->xyz;
}

int freesasa_coord_n(const freesasa_coord* c)
{
    assert(c);
    return c->n;
}

void freesasa_coord_translate(freesasa_coord *c, const double *xyz)
{
    assert(!c->is_const);
    assert(xyz);
    freesasa_coord_translate_xyz(c,xyz[0],xyz[1],xyz[2]);
}

void freesasa_coord_translate_xyz(freesasa_coord *c,
                                  double x, double y, double z)
{
    assert(c);
    assert(!c->is_const);

    for (int i = 0; i < c->n; ++i) {
        c->xyz[3*i]   += x;
        c->xyz[3*i+1] += y;
        c->xyz[3*i+2] += z;
    }
}

void freesasa_coord_scale(freesasa_coord *c, double s)
{
    assert(c);
    assert(!c->is_const);
    for (int i = 0; i < c->n*3; ++i) {
        c->xyz[i] *= s;
    }
}
