#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "freesasa_internal.h"
#include "nb.h"

#ifndef FREESASA_NB_CHUNK
#define FREESASA_NB_CHUNK 128
#endif

typedef struct cell cell;
struct cell {
    cell *nb[17]; /** includes self, only forward neighbors */
    int *atom;    /** indices of the atoms/coordinates in a cell */
    int n_nb;     /** number of neighbors to cell */
    int n_atoms;  /** number of atoms in cell */
};

static cell empty_cell = {{NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                           NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL},
                          NULL,
                          0,
                          0};

/** cell lists, divide space into boxes */
typedef struct cell_list {
    cell *cell;     /** the cells */
    int n;          /** number of cells */
    int nx, ny, nz; /** number of cells along each axis */
    double d;       /** cell size */
    double x_max, x_min;
    double y_max, y_min;
    double z_max, z_min;
} cell_list;

static struct cell_list empty_cell_list = {NULL, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

/** Finds the bounds of the cell list and writes them to the provided cell list */
static void
cell_list_bounds(cell_list *c,
                 const coord_t *coord)
{
    const int n = freesasa_coord_n(coord);
    int i;
    double d = c->d;
    const double *restrict v = freesasa_coord_i(coord, 0);
    double x = v[0], X = v[0], y = v[1], Y = v[1], z = v[2], Z = v[2];

    for (i = 1; i < n; ++i) {
        v = freesasa_coord_i(coord, i);
        x = fmin(v[0], x);
        X = fmax(v[0], X);
        y = fmin(v[1], y);
        Y = fmax(v[1], Y);
        z = fmin(v[2], z);
        Z = fmax(v[2], Z);
    }
    c->x_min = x - d / 2.;
    c->x_max = X + d / 2.;
    c->y_min = y - d / 2.;
    c->y_max = Y + d / 2.;
    c->z_min = z - d / 2.;
    c->z_max = Z + d / 2.;
    c->nx = (int)ceil((c->x_max - c->x_min) / d);
    c->ny = (int)ceil((c->y_max - c->y_min) / d);
    c->nz = (int)ceil((c->z_max - c->z_min) / d);
    c->n = c->nx * c->ny * c->nz;
}

static inline int
cell_index(const cell_list *c,
           int ix,
           int iy,
           int iz)
{
    assert(ix >= 0 && ix < c->nx);
    assert(iy >= 0 && iy < c->ny);
    return ix + c->nx * (iy + c->ny * iz);
}

/** Fill the neighbor list for a given cell, only "forward" neighbors considered */
static void
fill_nb(cell_list *c,
        int ix,
        int iy,
        int iz)
{
    cell *cell = &c->cell[cell_index(c, ix, iy, iz)];
    int n = 0, i, j, k;
    int xmin = ix > 0 ? ix - 1 : 0;
    int xmax = ix < c->nx - 1 ? ix + 1 : ix;
    int ymin = iy > 0 ? iy - 1 : 0;
    int ymax = iy < c->ny - 1 ? iy + 1 : iy;
    int zmin = iz > 0 ? iz - 1 : 0;
    int zmax = iz < c->nz - 1 ? iz + 1 : iz;
    for (i = xmin; i <= xmax; ++i) {
        for (j = ymin; j <= ymax; ++j) {
            for (k = zmin; k <= zmax; ++k) {
                /* Scalar product between (i-ix,j-iy,k-iz) and (1,1,1) should
                   be non-negative. Using only forward neighbors means
                   there's no double counting when comparing cells */
                if (i - ix + j - iy + k - iz >= 0) {
                    cell->nb[n] = &c->cell[cell_index(c, i, j, k)];
                    ++n;
                }
            }
        }
    }
    cell->n_nb = n;
    assert(n > 0);
}

/** find neighbors to all cells */
static void
get_nb(cell_list *c)
{
    int ix, iy, iz;

    for (ix = 0; ix < c->nx; ++ix) {
        for (iy = 0; iy < c->ny; ++iy) {
            for (iz = 0; iz < c->nz; ++iz) {
                fill_nb(c, ix, iy, iz);
            }
        }
    }
}

/** Get the cell index of a given atom */
static int
coord2cell_index(const cell_list *c,
                 const double *restrict xyz)
{
    double d = c->d;
    int ix = (int)((xyz[0] - c->x_min) / d);
    int iy = (int)((xyz[1] - c->y_min) / d);
    int iz = (int)((xyz[2] - c->z_min) / d);

    return cell_index(c, ix, iy, iz);
}

/**
   Assigns cells to each coordinate. Returns FREESASA_FAIL if realloc
   fails, FREESASA_SUCCESS else.
 */
static int
fill_cells(cell_list *c,
           const coord_t *coord)
{
    int i;
    cell *cell;
    int *a;
    const double *restrict v;

    for (i = 0; i < c->n; ++i) {
        c->cell[i].n_atoms = 0;
    }

    for (i = 0; i < freesasa_coord_n(coord); ++i) {
        v = freesasa_coord_i(coord, i);
        cell = &c->cell[coord2cell_index(c, v)];
        ++cell->n_atoms;
        a = cell->atom;
        cell->atom = realloc(cell->atom, sizeof(int) * cell->n_atoms);
        if (!cell->atom) {
            cell->atom = a;
            return mem_fail();
        }
        cell->atom[cell->n_atoms - 1] = i;
    }
    return FREESASA_SUCCESS;
}

/** Frees an object created by cell_list_new(). */
static void
cell_list_free(cell_list *c)
{
    int i;

    if (c) {
        if (c->cell) {
            for (i = 0; i < c->n; ++i)
                free(c->cell[i].atom);
        }
        free(c->cell);
        free(c);
    }
}

/**
    Creates a cell list with provided cell-size assigning cells to
    each of the provided coordinates. The created cell list should be
    freed using cell_list_free().

    Returns NULL if there are malloc fails.
 */
static cell_list *
cell_list_new(double cell_size,
              const coord_t *coord)
{
    int i;
    cell_list *c;

    assert(cell_size > 0);
    assert(coord);

    c = malloc(sizeof(cell_list));
    if (!c) {
        mem_fail();
        return NULL;
    }

    *c = empty_cell_list;

    c->d = cell_size;
    cell_list_bounds(c, coord);

    c->cell = malloc(sizeof(cell) * c->n);
    if (!c->cell) {
        cell_list_free(c);
        mem_fail();
        return NULL;
    }

    for (i = 0; i < c->n; ++i)
        c->cell[i] = empty_cell;

    if (fill_cells(c, coord)) {
        cell_list_free(c);
        mem_fail();
        return NULL;
    }

    get_nb(c);
    return c;
}

/** assumes max value in a is positive */
static double
max_array(const double *a,
          int n)
{
    int i;
    double max = 0;

    for (i = 0; i < n; ++i) {
        max = fmax(a[i], max);
    }

    return max;
}

/**
    Allocate memory for ::nb_list object. Tries to free everything
    and returns NULL if malloc somewhere along the way.
 */
static nb_list *
freesasa_nb_alloc(int n)
{
    int i;
    nb_list *nb;

    assert(n > 0);

    nb = malloc(sizeof(nb_list));
    if (!nb) {
        mem_fail();
        return NULL;
    }

    nb->n = n;

    /* in case the mallocs break, we can clean up in a safer way */
    nb->nn = NULL;
    nb->nb = NULL;
    nb->capacity = NULL;
    nb->xyd = nb->xd = nb->yd = NULL;

    nb->nn = malloc(sizeof(int) * n);
    nb->nb = malloc(sizeof(int *) * n);
    nb->xyd = malloc(sizeof(double *) * n);
    nb->xd = malloc(sizeof(double *) * n);
    nb->yd = malloc(sizeof(double *) * n);
    nb->capacity = malloc(sizeof(int) * n);

    if (!nb->nn || !nb->nb || !nb->xyd ||
        !nb->xd || !nb->yd || !nb->capacity) {
        free(nb->nn);
        free(nb->nb);
        free(nb->xyd);
        free(nb->xd);
        free(nb->yd);
        free(nb->capacity);
        free(nb);
        mem_fail();
        return NULL;
    }

    for (i = 0; i < n; ++i) {
        nb->nn[i] = 0;
        nb->capacity[i] = FREESASA_NB_CHUNK;
        /* again prepare for a potential cleanup */
        nb->nb[i] = NULL;
        nb->xyd[i] = nb->xd[i] = nb->yd[i] = NULL;
    }
    for (i = 0; i < n; ++i) {
        nb->nb[i] = malloc(sizeof(int) * FREESASA_NB_CHUNK);
        nb->xyd[i] = malloc(sizeof(double) * FREESASA_NB_CHUNK);
        nb->xd[i] = malloc(sizeof(double) * FREESASA_NB_CHUNK);
        nb->yd[i] = malloc(sizeof(double) * FREESASA_NB_CHUNK);
        if (!nb->nb[i] || !nb->xyd[i] || !nb->xd[i] || !nb->yd[i]) {
            freesasa_nb_free(nb);
            mem_fail();
            return NULL;
        }
    }
    return nb;
}

void freesasa_nb_free(nb_list *nb)
{
    int n, i;

    if (nb != NULL) {
        n = nb->n;
        if (nb->nb)
            for (i = 0; i < n; ++i)
                free(nb->nb[i]);
        if (nb->xyd)
            for (i = 0; i < n; ++i)
                free(nb->xyd[i]);
        if (nb->xd)
            for (i = 0; i < n; ++i)
                free(nb->xd[i]);
        if (nb->yd)
            for (i = 0; i < n; ++i)
                free(nb->yd[i]);
        free(nb->nb);
        free(nb->nn);
        free(nb->capacity);
        free(nb->xyd);
        free(nb->xd);
        free(nb->yd);
        free(nb);
    }
}

/**
    Increases sizes of arrays when they cross a threshold. Returns
    FREESASA_FAIL if realloc fails, FREESASA_SUCCESS else
 */
static int
chunk_up(nb_list *nb_list,
         int i)
{
    int nni = nb_list->nn[i];
    int **nbi, *nbi_b, new_cap;
    double **xydi, **xdi, **ydi, *xydi_b, *xdi_b, *ydi_b;

    if (nni > nb_list->capacity[i]) {
        nbi = &nb_list->nb[i];
        nbi_b = *nbi;
        xydi = &nb_list->xyd[i];
        xdi = &nb_list->xd[i];
        ydi = &nb_list->yd[i];
        xydi_b = *xydi;
        xdi_b = *xdi;
        ydi_b = *ydi;
        new_cap = (nb_list->capacity[i] += FREESASA_NB_CHUNK);

        *nbi = realloc(*nbi, sizeof(int) * new_cap);
        if (*nbi == NULL) {
            nb_list->nb[i] = nbi_b;
            return mem_fail();
        }

        *xydi = realloc(*xydi, sizeof(double) * new_cap);
        if (*xydi == NULL) {
            nb_list->xyd[i] = xydi_b;
            return mem_fail();
        }

        *xdi = realloc(*xdi, sizeof(double) * new_cap);
        if (*xdi == NULL) {
            nb_list->xd[i] = xdi_b;
            return mem_fail();
        }

        *ydi = realloc(*ydi, sizeof(double) * new_cap);
        if (*ydi == NULL) {
            nb_list->yd[i] = ydi_b;
            return mem_fail();
        }
    }
    return FREESASA_SUCCESS;
}

/**
    Assumes the coordinates i and j have been determined to be
    neighbors and adds them both to the provided nb lists,
    symmetrically.

    Returns FREESASA_FAIL if can't allocate memory. FREESASA_SUCCESS
    else.
*/
static int
nb_add_pair(nb_list *nb_list,
            int i,
            int j,
            double dx,
            double dy)
{
    int **nb;
    int *nn = nb_list->nn;
    int nni, nnj;
    double **xyd;
    double **xd;
    double **yd;
    double d;

    assert(i != j);

    nni = nn[i]++;
    nnj = nn[j]++;

    if (chunk_up(nb_list, i)) return mem_fail();
    if (chunk_up(nb_list, j)) return mem_fail();

    nb = nb_list->nb;
    xyd = nb_list->xyd;
    xd = nb_list->xd;
    yd = nb_list->yd;

    nb[i][nni] = j;
    nb[j][nnj] = i;

    d = sqrt(dx * dx + dy * dy);

    xyd[i][nni] = d;
    xyd[j][nnj] = d;

    xd[i][nni] = dx;
    xd[j][nnj] = -dx;
    yd[i][nni] = dy;
    yd[j][nnj] = -dy;

    return FREESASA_SUCCESS;
}

/**
    Fills the nb list for all contacts between coordinates
    belonging to the cells ci and cj. Handles the case ci == cj
    correctly.
*/
static int
nb_calc_cell_pair(nb_list *nb_list,
                  const coord_t *coord,
                  const double *radii,
                  const cell *ci,
                  const cell *cj)
{
    const double *restrict v = freesasa_coord_all(coord);
    double ri, rj, xi, yi, zi, xj, yj, zj,
        dx, dy, dz, cut2;
    int i, j, ia, ja;

    for (i = 0; i < ci->n_atoms; ++i) {
        ia = ci->atom[i];
        ri = radii[ia];
        xi = v[ia * 3];
        yi = v[ia * 3 + 1];
        zi = v[ia * 3 + 2];
        if (ci == cj)
            j = i + 1;
        else
            j = 0;
        /** the following loop is performance critical */
        for (; j < cj->n_atoms; ++j) {
            ja = cj->atom[j];
            rj = radii[ja];
            xj = v[ja * 3];
            yj = v[ja * 3 + 1];
            zj = v[ja * 3 + 2];
            cut2 = (ri + rj) * (ri + rj);
            dx = xj - xi;
            dy = yj - yi;
            dz = zj - zi;
            if (dx * dx + dy * dy + dz * dz < cut2) {
                if (nb_add_pair(nb_list, ia, ja, dx, dy))
                    return mem_fail();
            }
        }
    }
    return FREESASA_SUCCESS;
}

/**
    Iterates through the cells and records all contacts in the
    provided nb list
 */
static int
nb_fill_list(nb_list *nb_list,
             cell_list *c,
             const coord_t *coord,
             const double *radii)
{
    int nc = c->n, ic, jc;
    cell *ci, *cj;

    for (ic = 0; ic < nc; ++ic) {
        ci = &c->cell[ic];
        for (jc = 0; jc < ci->n_nb; ++jc) {
            cj = ci->nb[jc];
            if (nb_calc_cell_pair(nb_list, coord, radii, ci, cj))
                return mem_fail();
        }
    }
    return FREESASA_SUCCESS;
}

nb_list *
freesasa_nb_new(const coord_t *coord,
                const double *radii)
{
    double cell_size;
    cell_list *c;
    int n;
    nb_list *nb;

    if (coord == NULL || radii == NULL) return NULL;

    n = freesasa_coord_n(coord);
    nb = freesasa_nb_alloc(n);

    if (!nb) {
        mem_fail();
        return NULL;
    }

    cell_size = 2 * max_array(radii, n);
    assert(cell_size > 0);
    c = cell_list_new(cell_size, coord);
    if (c == NULL ||
        nb_fill_list(nb, c, coord, radii)) {
        mem_fail();
        freesasa_nb_free(nb);
        nb = NULL;
    }

    /* the cell lists are only a tool to generate the neighbor lists */
    cell_list_free(c);

    return nb;
}

int freesasa_nb_contact(const nb_list *nb,
                        int i,
                        int j)
{
    int k;
    assert(nb != NULL);
    assert(i < nb->n && i >= 0);
    assert(j < nb->n && j >= 0);

    for (k = 0; k < nb->nn[i]; ++k) {
        if (nb->nb[i][k] == j) return 1;
    }

    return 0;
}

#if USE_CHECK
#include <check.h>
#include <math.h>

START_TEST(test_cell)
{
    int na, i;
    const int n_atoms = 6;
    static const double v[] = {0, 0, 0, 1, 1, 1, -1, 1, -1, 2, 0, -2, 2, 2, 0, -5, 5, 5};
    static const double r[] = {4, 2, 2, 2, 2, 2};
    double r_max;
    cell_list *c;
    coord_t *coord = freesasa_coord_new();
    cell ci;

    freesasa_coord_append(coord, v, n_atoms);
    r_max = max_array(r, n_atoms);
    ck_assert(fabs(r_max - 4) < 1e-10);
    c = cell_list_new(r_max, coord);
    ck_assert(c != NULL);
    ck_assert(c->cell != NULL);
    ck_assert(fabs(c->d - r_max) < 1e-10);

    /* check bounds */
    ck_assert(c->x_min < -5);
    ck_assert(c->x_max > 2);
    ck_assert(c->y_min < 0);
    ck_assert(c->y_max > 5);
    ck_assert(c->z_min < -2);
    ck_assert(c->z_max > 5);

    /* check number of cells */
    ck_assert(c->nx * c->d >= 7);
    ck_assert(c->nx <= ceil(7 / r_max) + 1);
    ck_assert(c->ny * c->d >= 5);
    ck_assert(c->ny <= ceil(5 / r_max) + 1);
    ck_assert(c->nz * c->d >= 7);
    ck_assert(c->nz <= ceil(7 / r_max) + 1);
    ck_assert_int_eq(c->n, c->nx * c->ny * c->nz);

    /* check the individual cells */
    na = 0;
    ck_assert_int_eq(c->cell[0].n_nb, 8);
    ck_assert_int_eq(c->cell[c->n - 1].n_nb, 1);
    for (i = 0; i < c->n; ++i) {
        ci = c->cell[i];
        ck_assert(ci.n_atoms >= 0);
        if (ci.n_atoms > 0) ck_assert(ci.atom != NULL);
        ck_assert_int_ge(ci.n_nb, 1);
        ck_assert_int_le(ci.n_nb, 17);
        na += ci.n_atoms;
    }
    ck_assert_int_eq(na, n_atoms);
    cell_list_free(c);
    freesasa_coord_free(coord);
}
END_TEST

TCase *
test_nb_static()
{
    TCase *tc = tcase_create("nb.c static");
    tcase_add_test(tc, test_cell);

    return tc;
}

#endif /* USE_CHECK */
