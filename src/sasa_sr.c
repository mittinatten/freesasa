#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
#define _USE_MATH_DEFINES
#endif
#include <math.h>

#if USE_THREADS
#include <pthread.h>
#define MAX_SR_THREADS 16
#else
#define MAX_SR_THREADS 1
#endif

#include "freesasa_internal.h"
#include "nb.h"

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

/* calculation parameters (results stored in *sasa) */
typedef struct {
    int i1, i2; /* for multithreading, range of atoms */
    int thread_index;
    int n_atoms;
    int n_points;
    int n_threads;
    double probe_radius;
    const coord_t *xyz;
    coord_t *srp;                      /* test-points */
    coord_t *tp_local[MAX_SR_THREADS]; /* coord object for storing intermediates */
    int *spcount[MAX_SR_THREADS];
    double *r;
    double *r2;
    nb_list *nb;
    double *sasa;
} sr_data;

#if USE_THREADS
static int sr_do_threads(int n_threads, sr_data *sr);
static void *sr_thread(void *arg);
#endif

static double
sr_atom_area(int i, const sr_data *sr, int thread_index) __attrib_pure__;

static coord_t *
test_points(int N)
{
    /* Golden section spiral on a sphere
       from http://web.archive.org/web/20120421191837/http://www.cgafaq.info/wiki/Evenly_distributed_points_on_sphere */
    double dlong = M_PI * (3 - sqrt(5)), dz = 2.0 / N, longitude = 0, z = 1 - dz / 2, r;
    coord_t *coord = freesasa_coord_new();
    double *tp = malloc(3 * N * sizeof(double)), *p;
    if (tp == NULL || coord == NULL) {
        mem_fail();
        goto cleanup;
    }

    for (p = tp; p - tp < 3 * N; p += 3) {
        r = sqrt(1 - z * z);
        p[0] = cos(longitude) * r;
        p[1] = sin(longitude) * r;
        p[2] = z;
        z -= dz;
        longitude += dlong;
    }

    if (freesasa_coord_append(coord, tp, N) == FREESASA_FAIL) {
        fail_msg("");
        goto cleanup;
    }
    free(tp);

    return coord;

cleanup:
    free(tp);
    freesasa_coord_free(coord);
    return NULL;
}

/* free contents */
void release_sr(sr_data *sr)
{
    int i;

    freesasa_coord_free(sr->srp);
    freesasa_nb_free(sr->nb);
    free(sr->r);
    free(sr->r2);

    for (i = 0; i < sr->n_threads; ++i) {
        freesasa_coord_free(sr->tp_local[i]);
        free(sr->spcount[i]);
    }
}

int init_sr(sr_data *sr,
            double *sasa,
            const coord_t *xyz,
            const double *r,
            double probe_radius,
            int n_points,
            int n_threads)
{
    int n_atoms = freesasa_coord_n(xyz), i;
    coord_t *srp = test_points(n_points);
    double ri;

    if (srp == NULL) return fail_msg("failed to initialize test points");

    /* store parameters and reference arrays */
    sr->n_atoms = n_atoms;
    sr->n_points = n_points;
    sr->n_threads = n_threads;
    sr->probe_radius = probe_radius;
    sr->xyz = xyz;
    sr->srp = srp;
    sr->sasa = sasa;
    sr->nb = NULL;

    /* should be done before any mallocs (to avoid problems in potential cleanup) */
    for (i = 0; i < n_threads; ++i) {
        sr->tp_local[i] = NULL;
        sr->spcount[i] = NULL;
    }

    sr->r = malloc(sizeof(double) * n_atoms);
    sr->r2 = malloc(sizeof(double) * n_atoms);

    if (sr->r == NULL || sr->r2 == NULL) goto cleanup;

    for (i = 0; i < n_atoms; ++i) {
        ri = r[i] + probe_radius;
        sr->r[i] = ri;
        sr->r2[i] = ri * ri;
    }

    for (i = 0; i < n_threads; ++i) {
        sr->tp_local[i] = freesasa_coord_clone(sr->srp);
        sr->spcount[i] = malloc(sizeof(int) * n_points);
        if (sr->tp_local[i] == NULL || sr->spcount[i] == NULL) {
            goto cleanup;
        }
    }

    /* calculate distances */
    sr->nb = freesasa_nb_new(xyz, sr->r);
    if (sr->nb == NULL) goto cleanup;

    return FREESASA_SUCCESS;

cleanup:
    release_sr(sr);
    return mem_fail();
}

int freesasa_shrake_rupley(double *sasa,
                           const coord_t *xyz,
                           const double *r,
                           const freesasa_parameters *param)
{
    int i, n_atoms, n_threads = param->n_threads, resolution, return_value;
    double probe_radius = param->probe_radius;
    sr_data sr;

    assert(sasa);
    assert(xyz);
    assert(r);

    if (param == NULL) param = &freesasa_default_parameters;

    n_atoms = freesasa_coord_n(xyz);
    n_threads = param->n_threads;
    resolution = param->shrake_rupley_n_points;
    return_value = FREESASA_SUCCESS;

    if (n_threads > MAX_SR_THREADS) {
        return fail_msg("S&R does not support more than %d threads", MAX_SR_THREADS);
    }
    if (resolution <= 0) {
        return fail_msg("%f test points invalid resolution in S&R, must be > 0\n", resolution);
    }
    if (n_atoms == 0) return freesasa_warn("in %s(): empty coordinates", __func__);
    if (n_threads > n_atoms) {
        n_threads = n_atoms;
        freesasa_warn("no sense in having more threads than atoms, only using %d threads",
                      n_threads);
    }

    if (init_sr(&sr, sasa, xyz, r, probe_radius, resolution, n_threads))
        return FREESASA_FAIL;

    /* calculate SASA */
    if (n_threads > 1) {
#if USE_THREADS
        return_value = sr_do_threads(n_threads, &sr);
#else
        return_value = freesasa_warn("in %s(): program compiled for single-threaded use, "
                                     "but multiple threads were requested, will "
                                     "proceed in single-threaded mode\n",
                                     __func__);
        n_threads = 1;
#endif
    }
    if (n_threads == 1) {
        /* don't want the overhead of generating threads if only one is used */
        for (i = 0; i < n_atoms; ++i) {
            sasa[i] = sr_atom_area(i, &sr, 0);
        }
    }
    release_sr(&sr);
    return return_value;
}

#if USE_THREADS
static int
sr_do_threads(int n_threads,
              sr_data *sr)
{
    pthread_t thread[MAX_SR_THREADS];
    sr_data srt[MAX_SR_THREADS];
    int thread_block_size = sr->n_atoms / n_threads;
    int res, return_value = FREESASA_SUCCESS;
    int threads_created = 0, t;

    /*  divide atoms evenly over threads */
    for (t = 0; t < n_threads; ++t) {
        srt[t] = *sr;
        srt[t].i1 = t * thread_block_size;
        if (t == n_threads - 1)
            srt[t].i2 = sr->n_atoms;
        else
            srt[t].i2 = (t + 1) * thread_block_size;
        srt[t].thread_index = t;
        res = pthread_create(&thread[t], NULL, sr_thread, (void *)&srt[t]);
        if (res) {
            return_value = fail_msg(freesasa_thread_error(res));
            break;
        }
        ++threads_created;
    }
    for (t = 0; t < threads_created; ++t) {
        res = pthread_join(thread[t], NULL);
        if (res) {
            return_value = fail_msg(freesasa_thread_error(res));
        }
    }
    return return_value;
}

static void *
sr_thread(void *arg)
{
    int i;
    sr_data *sr = ((sr_data *)arg);

    for (i = sr->i1; i < sr->i2; ++i) {
        /* mutex should not be necessary, writes to non-overlapping regions */
        sr->sasa[i] = sr_atom_area(i, sr, sr->thread_index);
    }
    pthread_exit(NULL);
}
#endif

static double
sr_atom_area(int i,
             const sr_data *sr,
             int thread_index)
{
    const int n_points = sr->n_points;
    /* this array keeps track of which testpoints belonging to
       a certain atom do not overlap with any other atoms */
    int *spcount = sr->spcount[thread_index];
    const int nni = sr->nb->nn[i];
    const int *restrict nbi = sr->nb->nb[i];
    const double ri = sr->r[i];
    const double *restrict r2 = sr->r2;
    const double *restrict v = freesasa_coord_all(sr->xyz);
    const double *restrict vi = v + 3 * i;
    const double *restrict tp;
    int n_surface = 0, current_nb, a, j, k;
    double dx, dy, dz;
    /* testpoints for this atom */
    coord_t *restrict tp_coord_ri = sr->tp_local[thread_index];

    freesasa_coord_copy(tp_coord_ri, sr->srp);
    freesasa_coord_scale(tp_coord_ri, ri);
    freesasa_coord_translate(tp_coord_ri, vi);
    tp = freesasa_coord_all(tp_coord_ri);

    /* initialize with all surface points hidden */
    memset(spcount, 0, n_points * sizeof(int));

    /* Using the trick from NSOL to check points one by one for all
       atoms, start comparing with the first neighbor. If there is no
       overlap for a given test-point, try with other neighbors
       instead. Would probably work even better if test points were
       organized in patches and not spirals. */
    current_nb = 0;
    for (j = 0; j < n_points; ++j) {
        /* a is the index of the atom under consideration */
        a = nbi[current_nb];
        dx = tp[j * 3] - v[a * 3];
        dy = tp[j * 3 + 1] - v[a * 3 + 1];
        dz = tp[j * 3 + 2] - v[a * 3 + 2];
        if (dx * dx + dy * dy + dz * dz > r2[a]) {
            k = 0;
            for (; k < nni; ++k) {
                a = nbi[k];
                dx = tp[j * 3] - v[a * 3];
                dy = tp[j * 3 + 1] - v[a * 3 + 1];
                dz = tp[j * 3 + 2] - v[a * 3 + 2];
                if (dx * dx + dy * dy + dz * dz <= r2[a]) {
                    current_nb = k;
                    break;
                }
            }
            /* we have gone through the whole list without overlap */
            if (k == nni) spcount[j] = 1;
        }
    }
    for (k = 0; k < n_points; ++k) {
        if (spcount[k]) ++n_surface;
    }

    return (4.0 * M_PI * ri * ri * n_surface) / n_points;
}
