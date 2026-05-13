#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#if USE_OPENMP
#include <omp.h>
#endif

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

/* ─────────────────────────────────────────────────────────────────────
   Thread-local pair buffer for parallel neighbor list construction
   ───────────────────────────────────────────────────────────────────── */

/** A neighbor pair discovered during cell scanning */
typedef struct {
    int i, j;      /** atom indices */
    double dx, dy; /** signed displacements */
} nb_pair;

/** Growable buffer of pairs for one thread */
typedef struct {
    nb_pair *pairs;
    int n;
    int capacity;
} pair_buffer;

static int
pair_buffer_init(pair_buffer *pb, int initial_cap)
{
    pb->pairs = malloc(sizeof(nb_pair) * initial_cap);
    if (!pb->pairs) return mem_fail();
    pb->n = 0;
    pb->capacity = initial_cap;
    return FREESASA_SUCCESS;
}

static void
pair_buffer_free(pair_buffer *pb)
{
    free(pb->pairs);
    pb->pairs = NULL;
    pb->n = 0;
    pb->capacity = 0;
}

static int
pair_buffer_add(pair_buffer *pb, int i, int j, double dx, double dy)
{
    if (pb->n >= pb->capacity) {
        int new_cap = pb->capacity * 2;
        nb_pair *tmp = realloc(pb->pairs, sizeof(nb_pair) * new_cap);
        if (!tmp) return mem_fail();
        pb->pairs = tmp;
        pb->capacity = new_cap;
    }
    nb_pair *p = &pb->pairs[pb->n++];
    p->i = i;
    p->j = j;
    p->dx = dx;
    p->dy = dy;
    return FREESASA_SUCCESS;
}

/* ─────────────────────────────────────────────────────────────────────
   Parallel neighbor list construction
   ───────────────────────────────────────────────────────────────────── */

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

    nb->nn = calloc(n, sizeof(int));
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
        nb->capacity[i] = 0;
        /* prepare for a potential cleanup */
        nb->nb[i] = NULL;
        nb->xyd[i] = nb->xd[i] = nb->yd[i] = NULL;
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
    Scan cell pairs for a range of cells and collect neighbor pairs
    into a thread-local buffer. No writes to shared data.
*/
static int
nb_scan_cells_range(pair_buffer *pb,
                    cell_list *c,
                    const coord_t *coord,
                    const double *radii,
                    int cell_start,
                    int cell_end)
{
    const double *restrict v = freesasa_coord_all(coord);
    int ic, jc;
    double ri, rj, xi, yi, zi, xj, yj, zj, dx, dy, dz, cut2;
    int i, j, ia, ja;
    cell *ci, *cj;

    for (ic = cell_start; ic < cell_end; ++ic) {
        ci = &c->cell[ic];
        for (jc = 0; jc < ci->n_nb; ++jc) {
            cj = ci->nb[jc];
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
                        if (pair_buffer_add(pb, ia, ja, dx, dy))
                            return mem_fail();
                    }
                }
            }
        }
    }
    return FREESASA_SUCCESS;
}

/**
    Build neighbor list from pair buffers.

    Phase 1: Count neighbors for each atom (parallel-safe since each
             pair contributes to two atoms)
    Phase 2: Allocate per-atom arrays
    Phase 3: Fill per-atom arrays from pairs
*/
static int
nb_build_from_pairs(nb_list *nb,
                    pair_buffer *buffers,
                    int n_buffers)
{
    int b, p, n = nb->n;
    int *nn = nb->nn;

    /* Phase 1: count neighbors per atom */
    for (b = 0; b < n_buffers; ++b) {
        pair_buffer *pb = &buffers[b];
        for (p = 0; p < pb->n; ++p) {
            nn[pb->pairs[p].i]++;
            nn[pb->pairs[p].j]++;
        }
    }

    /* Phase 2: allocate per-atom arrays */
    {
        int i;
        for (i = 0; i < n; ++i) {
            int cap = nn[i] > 0 ? nn[i] : 1;
            nb->capacity[i] = cap;
            nb->nb[i] = malloc(sizeof(int) * cap);
            nb->xyd[i] = malloc(sizeof(double) * cap);
            nb->xd[i] = malloc(sizeof(double) * cap);
            nb->yd[i] = malloc(sizeof(double) * cap);
            if (!nb->nb[i] || !nb->xyd[i] || !nb->xd[i] || !nb->yd[i]) {
                return mem_fail();
            }
            nn[i] = 0; /* reset counts for fill phase */
        }
    }

    /* Phase 3: fill per-atom arrays from pairs (symmetric) */
    for (b = 0; b < n_buffers; ++b) {
        pair_buffer *pb = &buffers[b];
        for (p = 0; p < pb->n; ++p) {
            int i = pb->pairs[p].i;
            int j = pb->pairs[p].j;
            double dx = pb->pairs[p].dx;
            double dy = pb->pairs[p].dy;
            double d = sqrt(dx * dx + dy * dy);
            int nni = nn[i]++;
            int nnj = nn[j]++;

            nb->nb[i][nni] = j;
            nb->nb[j][nnj] = i;

            nb->xyd[i][nni] = d;
            nb->xyd[j][nnj] = d;

            nb->xd[i][nni] = dx;
            nb->xd[j][nnj] = -dx;

            nb->yd[i][nni] = dy;
            nb->yd[j][nnj] = -dy;
        }
    }

    return FREESASA_SUCCESS;
}

/**
    Parallelized neighbor list fill using OpenMP.
    Falls back to serial if OpenMP is not available or n_cells is small.
*/
static int
nb_fill_list_parallel(nb_list *nb,
                      cell_list *c,
                      const coord_t *coord,
                      const double *radii)
{
    int nc = c->n;
    int ret = FREESASA_SUCCESS;

#if USE_OPENMP
    int n_threads = omp_get_max_threads();
    if (n_threads < 1) n_threads = 1;
    /* For very small cell lists, don't bother with parallelism overhead */
    if (nc < 64) n_threads = 1;
#else
    int n_threads = 1;
#endif

    pair_buffer *buffers = calloc(n_threads, sizeof(pair_buffer));
    if (!buffers) return mem_fail();

    {
        int t;
        for (t = 0; t < n_threads; ++t) {
            if (pair_buffer_init(&buffers[t], 4096)) {
                int u;
                for (u = 0; u < t; ++u) pair_buffer_free(&buffers[u]);
                free(buffers);
                return mem_fail();
            }
        }
    }

#if USE_OPENMP
    if (n_threads > 1) {
        int shared_error = 0;
        #pragma omp parallel num_threads(n_threads) default(none) \
            shared(buffers, c, coord, radii, nc, shared_error)
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
            int cells_per_thread = (nc + nthreads - 1) / nthreads;
            int start = tid * cells_per_thread;
            int end = start + cells_per_thread;
            if (end > nc) end = nc;
            if (start < nc) {
                if (nb_scan_cells_range(&buffers[tid], c, coord, radii,
                                        start, end)) {
                    #pragma omp atomic write
                    shared_error = 1;
                }
            }
        }
        if (shared_error) {
            ret = FREESASA_FAIL;
            goto cleanup;
        }
    } else
#endif
    {
        /* Single-threaded path */
        if (nb_scan_cells_range(&buffers[0], c, coord, radii, 0, nc)) {
            ret = FREESASA_FAIL;
            goto cleanup;
        }
    }

    /* Merge phase: build the neighbor list from collected pairs */
    if (nb_build_from_pairs(nb, buffers, n_threads)) {
        ret = FREESASA_FAIL;
    }

cleanup:
    {
        int t;
        for (t = 0; t < n_threads; ++t) {
            pair_buffer_free(&buffers[t]);
        }
        free(buffers);
    }
    return ret;
}

/* Legacy serial fill for reference and fallback */
static int
nb_add_pair_serial(nb_list *nb_list,
                   int i,
                   int j,
                   double dx,
                   double dy)
{
    int *nn = nb_list->nn;
    int nni, nnj;
    double d;

    assert(i != j);

    nni = nn[i]++;
    nnj = nn[j]++;

    /* grow arrays if needed */
    if (nni >= nb_list->capacity[i]) {
        int new_cap = nb_list->capacity[i] + FREESASA_NB_CHUNK;
        nb_list->capacity[i] = new_cap;
        nb_list->nb[i] = realloc(nb_list->nb[i], sizeof(int) * new_cap);
        nb_list->xyd[i] = realloc(nb_list->xyd[i], sizeof(double) * new_cap);
        nb_list->xd[i] = realloc(nb_list->xd[i], sizeof(double) * new_cap);
        nb_list->yd[i] = realloc(nb_list->yd[i], sizeof(double) * new_cap);
        if (!nb_list->nb[i] || !nb_list->xyd[i] || !nb_list->xd[i] || !nb_list->yd[i])
            return mem_fail();
    }
    if (nnj >= nb_list->capacity[j]) {
        int new_cap = nb_list->capacity[j] + FREESASA_NB_CHUNK;
        nb_list->capacity[j] = new_cap;
        nb_list->nb[j] = realloc(nb_list->nb[j], sizeof(int) * new_cap);
        nb_list->xyd[j] = realloc(nb_list->xyd[j], sizeof(double) * new_cap);
        nb_list->xd[j] = realloc(nb_list->xd[j], sizeof(double) * new_cap);
        nb_list->yd[j] = realloc(nb_list->yd[j], sizeof(double) * new_cap);
        if (!nb_list->nb[j] || !nb_list->xyd[j] || !nb_list->xd[j] || !nb_list->yd[j])
            return mem_fail();
    }

    nb_list->nb[i][nni] = j;
    nb_list->nb[j][nnj] = i;

    d = sqrt(dx * dx + dy * dy);

    nb_list->xyd[i][nni] = d;
    nb_list->xyd[j][nnj] = d;

    nb_list->xd[i][nni] = dx;
    nb_list->xd[j][nnj] = -dx;
    nb_list->yd[i][nni] = dy;
    nb_list->yd[j][nnj] = -dy;

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
        nb_fill_list_parallel(nb, c, coord, radii)) {
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
