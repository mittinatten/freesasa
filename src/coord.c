#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "coord.h"
#include "freesasa_internal.h"

coord_t *
freesasa_coord_new()
{
    coord_t *c = malloc(sizeof(coord_t));
    if (c != NULL) {
        c->xyz = NULL;
        c->n = 0;
        c->is_linked = 0;
    } else {
        mem_fail();
    }
    return c;
}

void freesasa_coord_free(coord_t *c)
{
    if (c) {
        if (c->xyz && !c->is_linked) free(c->xyz);
        free(c);
    }
}
static void
coord_clear(coord_t *c)
{
    assert(c);
    assert(!c->is_linked);
    if (c->xyz) {
        free(c->xyz);
        c->xyz = NULL;
    }
    c->n = 0;
}

coord_t *
freesasa_coord_clone(const coord_t *src)
{
    coord_t *c = freesasa_coord_new();

    assert(src);

    if (c == NULL) {
        mem_fail();
        return NULL;
    }

    if (freesasa_coord_set_all(c, src->xyz, src->n) != FREESASA_SUCCESS) {
        fail_msg("");
        return NULL;
    }

    return c;
}

int freesasa_coord_copy(coord_t *target, const coord_t *src)
{
    if (src->n != target->n) return FREESASA_FAIL;
    memcpy(target->xyz, src->xyz, sizeof(double) * src->n * 3);
    return FREESASA_SUCCESS;
}

coord_t *
freesasa_coord_new_linked(const double *xyz,
                          int n)
{
    coord_t *c = freesasa_coord_new();

    assert(xyz);
    assert(n > 0);

    if (c != NULL) {
        c->xyz = (double *)xyz;
        c->n = n;
        c->is_linked = 1;
    } else
        mem_fail();
    return c;
}

int freesasa_coord_append(coord_t *c,
                          const double *xyz,
                          int n)
{
    int n_old;
    double *xyz_old;

    assert(c);
    assert(xyz);
    assert(!c->is_linked);

    n_old = c->n;
    xyz_old = c->xyz;

    if (n == 0) return FREESASA_SUCCESS;

    c->xyz = (double *)realloc(c->xyz, sizeof(double) * 3 * (c->n + n));
    if (c->xyz == NULL) {
        free(xyz_old);
        return mem_fail();
    }
    c->n += n;

    memcpy(&(c->xyz[3 * n_old]), xyz, sizeof(double) * n * 3);
    return FREESASA_SUCCESS;
}

int freesasa_coord_append_xyz(coord_t *c,
                              const double *x,
                              const double *y,
                              const double *z, int n)
{
    double *xyz;
    int i;

    assert(c);
    assert(x);
    assert(y);
    assert(z);
    assert(!c->is_linked);

    if (n == 0) return FREESASA_SUCCESS;

    xyz = malloc(sizeof(double) * n * 3);
    if (xyz == NULL) return mem_fail();

    for (i = 0; i < n; ++i) {
        xyz[i * 3] = x[i];
        xyz[i * 3 + 1] = y[i];
        xyz[i * 3 + 2] = z[i];
    }
    if (freesasa_coord_append(c, xyz, n) == FREESASA_SUCCESS) {
        free(xyz);
        return FREESASA_SUCCESS;
    }
    return mem_fail();
}

void freesasa_coord_set_i(coord_t *c,
                          int i,
                          const double *xyz)
{
    assert(c);
    assert(xyz);
    assert(i < c->n && i >= 0);
    assert(!c->is_linked);

    memcpy(&c->xyz[i * 3], xyz, 3 * sizeof(double));
}

void freesasa_coord_set_i_xyz(coord_t *c,
                              int i,
                              double x,
                              double y,
                              double z)
{
    double *v_i;

    assert(c);
    assert(c->n > i);
    assert(i >= 0);
    assert(!c->is_linked);

    v_i = &c->xyz[i * 3];

    *(v_i++) = x;
    *(v_i++) = y;
    *v_i = z;
}

int freesasa_coord_set_all(coord_t *c,
                           const double *xyz,
                           int n)
{
    int result;

    assert(c);
    assert(xyz);

    coord_clear(c);
    result = freesasa_coord_append(c, xyz, n);

    if (result != FREESASA_SUCCESS) {
        fail_msg("");
    }

    return result;
}

int freesasa_coord_set_all_xyz(coord_t *c,
                               const double *x,
                               const double *y,
                               const double *z,
                               int n)
{
    assert(c);
    assert(x);
    assert(y);
    assert(z);
    coord_clear(c);
    return freesasa_coord_append_xyz(c, x, y, z, n);
}

void freesasa_coord_set_length_i(coord_t *c,
                                 int i,
                                 double l)
{
    double x, y, z, r;

    assert(c);
    assert(c->xyz);
    assert(!c->is_linked);
    assert(i >= 0 && i < c->n);
    assert(l >= 0);

    x = c->xyz[3 * i];
    y = c->xyz[3 * i + 1];
    z = c->xyz[3 * i + 2];
    r = sqrt(x * x + y * y + z * z);
    c->xyz[3 * i] *= l / r;
    c->xyz[3 * i + 1] *= l / r;
    c->xyz[3 * i + 2] *= l / r;
}

void freesasa_coord_set_length_all(coord_t *c,
                                   double l)
{
    int i;

    assert(c);
    assert(!c->is_linked);

    for (i = 0; i < c->n; ++i)
        freesasa_coord_set_length_i(c, i, l);
}

const double *
freesasa_coord_i(const coord_t *c,
                 int i)
{
    assert(c);
    assert(i < c->n);
    assert(i >= 0);
    return &c->xyz[3 * i];
}

double
freesasa_coord_dist(const coord_t *c,
                    int i,
                    int j)
{
    return sqrt(freesasa_coord_dist2(c, i, j));
}

static inline double
dist2(const double *v1,
      const double *v2)
{
    double dx = v1[0] - v2[0], dy = v1[1] - v2[1], dz = v1[2] - v2[2];
    return dx * dx + dy * dy + dz * dz;
}

double
freesasa_coord_dist2(const coord_t *c,
                     int i,
                     int j)
{
    double *v1 = &c->xyz[3 * i];
    double *v2 = &c->xyz[3 * j];
    return dist2(v1, v2);
}

double
freesasa_coord_dist2_12(const coord_t *c1,
                        const coord_t *c2,
                        int i1,
                        int i2)
{
    double *v1 = &c1->xyz[3 * i1];
    double *v2 = &c2->xyz[3 * i2];
    return dist2(v1, v2);
}

const double *
freesasa_coord_all(const coord_t *c)
{
    assert(c);
    return c->xyz;
}

int freesasa_coord_n(const coord_t *c)
{
    assert(c);
    return c->n;
}

void freesasa_coord_translate(coord_t *c,
                              const double *xyz)
{
    assert(!c->is_linked);
    assert(xyz);
    freesasa_coord_translate_xyz(c, xyz[0], xyz[1], xyz[2]);
}

void freesasa_coord_translate_xyz(coord_t *c,
                                  double x,
                                  double y,
                                  double z)
{
    int i;

    assert(c);
    assert(!c->is_linked);

    for (i = 0; i < c->n; ++i) {
        c->xyz[3 * i] += x;
        c->xyz[3 * i + 1] += y;
        c->xyz[3 * i + 2] += z;
    }
}

void freesasa_coord_scale(coord_t *c,
                          double s)
{
    int i;

    assert(c);
    assert(!c->is_linked);

    for (i = 0; i < c->n * 3; ++i) {
        c->xyz[i] *= s;
    }
}
