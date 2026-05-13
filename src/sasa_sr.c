/* sasa_sr.c — Shrake-Rupley SASA with OpenMP parallelism and SIMD vectorization.
 *
 * Parallelism layers:
 *   1. Atom-level: #pragma omp parallel for distributes atoms across threads.
 *   2. SIMD inner loop: neighbor distance checks vectorized with AVX-512 floats
 *      (16-wide), AVX2 (8-wide), or auto-vectorization fallback.
 *
 * Key data-layout change: per-atom SoA neighbor cache (float32)
 *   nb_cx[k], nb_cy[k], nb_cz[k], nb_r2[k] for k=0..nni-1
 * This gives 2x SIMD density vs double and ensures sequential memory access.
 */

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

#if USE_OPENMP
#include <omp.h>
#endif

/* SIMD intrinsics — only include on x86 with AVX2/AVX-512 */
#if defined(__AVX512F__) && defined(__AVX512DQ__)
#  include <immintrin.h>
#  define FREESASA_USE_AVX512 1
#elif defined(__AVX2__)
#  include <immintrin.h>
#  define FREESASA_USE_AVX2 1
#endif

#include "freesasa_internal.h"
#include "nb.h"

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

/*
 * Portable aligned memory allocation (C99-compatible).
 *
 * aligned_alloc() is C11 only.  On POSIX systems (Linux, macOS) we use
 * posix_memalign(); on MSVC we use _aligned_malloc()/_aligned_free().
 * The returned pointer must be released with freesasa_aligned_free().
 */
#if defined(_MSC_VER)
#  include <malloc.h>
static void *freesasa_aligned_malloc(size_t alignment, size_t size)
{
    return _aligned_malloc(size, alignment);
}
static void freesasa_aligned_free(void *p) { _aligned_free(p); }
#else
#  include <stdlib.h>
static void *freesasa_aligned_malloc(size_t alignment, size_t size)
{
    void *p = NULL;
    if (posix_memalign(&p, alignment, size) != 0) return NULL;
    return p;
}
static void freesasa_aligned_free(void *p) { free(p); }
#endif

/* Per-atom neighbor SoA cache (float32 for SIMD width).
 * For each atom i, arrays of length nb->nn[i] are allocated in a
 * single 64-byte-aligned block so SIMD loads can use vmovaps directly.
 */
typedef struct {
    float *x;   /* neighbor x-coords */
    float *y;
    float *z;
    float *r2;  /* squared (radius + probe)  */
    int    n;   /* number of neighbors       */
} nb_cache_t;

/* calculation parameters (results stored in *sasa) */
typedef struct {
    int n_atoms;
    int n_points;
    int n_threads;
    double probe_radius;
    const coord_t *xyz;
    coord_t *srp;         /* unit test-points */
    coord_t **tp_local;   /* per-thread scaled+translated test points */
    int **spcount;         /* per-thread hidden-point flags */
    double *r;
    double *r2;
    nb_list *nb;
    nb_cache_t *nb_cache; /* per-atom SoA neighbor cache [n_atoms] */
    double *sasa;
} sr_data;

static double
sr_atom_area(int i, const sr_data *sr, int thread_index) __attrib_pure__;

/* ──────────────────────────────────────────────────────────────────────────
 * Golden-section spiral test points
 * ────────────────────────────────────────────────────────────────────────── */
static coord_t *
test_points(int N)
{
    double dlong = M_PI * (3 - sqrt(5)), dz = 2.0 / N, longitude = 0,
           z = 1 - dz / 2, r;
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

/* ──────────────────────────────────────────────────────────────────────────
 * Build per-atom SoA neighbor cache in float32
 * ────────────────────────────────────────────────────────────────────────── */
static int
sr_build_nb_cache(sr_data *sr)
{
    int n_atoms = sr->n_atoms;
    const double *v = freesasa_coord_all(sr->xyz);
    const double *r2 = sr->r2;
    const nb_list *nb = sr->nb;

    sr->nb_cache = malloc(sizeof(nb_cache_t) * n_atoms);
    if (!sr->nb_cache) return mem_fail();

    for (int i = 0; i < n_atoms; ++i) {
        int nni = nb->nn[i];
        const int *nbi = nb->nb[i];
        nb_cache_t *c = &sr->nb_cache[i];
        c->n = nni;

        if (nni == 0) {
            c->x = c->y = c->z = c->r2 = NULL;
            continue;
        }

        /* Pad to nearest multiple of 16 for AVX-512 alignment */
        int padded = (nni + 15) & ~15;
        float *buf = (float *)freesasa_aligned_malloc(64, sizeof(float) * padded * 4);
        if (!buf) { free(sr->nb_cache); sr->nb_cache = NULL; return mem_fail(); }

        c->x  = buf;
        c->y  = buf + padded;
        c->z  = buf + padded * 2;
        c->r2 = buf + padded * 3;

        for (int k = 0; k < nni; ++k) {
            int a = nbi[k];
            c->x[k]  = (float)v[3 * a];
            c->y[k]  = (float)v[3 * a + 1];
            c->z[k]  = (float)v[3 * a + 2];
            c->r2[k] = (float)r2[a];
        }
        /* Pad tail with impossible values so distance check always fails */
        for (int k = nni; k < padded; ++k) {
            c->x[k] = c->y[k] = c->z[k] = 1e18f;
            c->r2[k] = 0.0f;
        }
    }
    return FREESASA_SUCCESS;
}

static void
sr_free_nb_cache(sr_data *sr)
{
    if (!sr->nb_cache) return;
    for (int i = 0; i < sr->n_atoms; ++i) {
        if (sr->nb_cache[i].x) freesasa_aligned_free(sr->nb_cache[i].x);
    }
    free(sr->nb_cache);
    sr->nb_cache = NULL;
}

/* ──────────────────────────────────────────────────────────────────────────
 * Release sr_data
 * ────────────────────────────────────────────────────────────────────────── */
void release_sr(sr_data *sr)
{
    int i;
    freesasa_coord_free(sr->srp);
    freesasa_nb_free(sr->nb);
    free(sr->r);
    free(sr->r2);
    sr_free_nb_cache(sr);

    if (sr->tp_local) {
        for (i = 0; i < sr->n_threads; ++i)
            freesasa_coord_free(sr->tp_local[i]);
        free(sr->tp_local);
    }
    if (sr->spcount) {
        for (i = 0; i < sr->n_threads; ++i)
            free(sr->spcount[i]);
        free(sr->spcount);
    }
}

/* ──────────────────────────────────────────────────────────────────────────
 * Initialize sr_data
 * ────────────────────────────────────────────────────────────────────────── */
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

    sr->n_atoms     = n_atoms;
    sr->n_points    = n_points;
    sr->n_threads   = n_threads;
    sr->probe_radius = probe_radius;
    sr->xyz         = xyz;
    sr->srp         = srp;
    sr->sasa        = sasa;
    sr->nb          = NULL;
    sr->tp_local    = NULL;
    sr->spcount     = NULL;
    sr->nb_cache    = NULL;

    sr->tp_local = calloc(n_threads, sizeof(coord_t *));
    sr->spcount  = calloc(n_threads, sizeof(int *));
    if (!sr->tp_local || !sr->spcount) goto cleanup;

    sr->r  = malloc(sizeof(double) * n_atoms);
    sr->r2 = malloc(sizeof(double) * n_atoms);
    if (!sr->r || !sr->r2) goto cleanup;

    for (i = 0; i < n_atoms; ++i) {
        ri = r[i] + probe_radius;
        sr->r[i]  = ri;
        sr->r2[i] = ri * ri;
    }

    for (i = 0; i < n_threads; ++i) {
        sr->tp_local[i] = freesasa_coord_clone(sr->srp);
        sr->spcount[i]  = malloc(sizeof(int) * n_points);
        if (!sr->tp_local[i] || !sr->spcount[i]) goto cleanup;
    }

    sr->nb = freesasa_nb_new(xyz, sr->r);
    if (!sr->nb) goto cleanup;

    /* Build per-atom float32 SoA neighbor cache for SIMD */
    if (sr_build_nb_cache(sr) != FREESASA_SUCCESS) goto cleanup;

    return FREESASA_SUCCESS;
cleanup:
    release_sr(sr);
    return mem_fail();
}

/* ──────────────────────────────────────────────────────────────────────────
 * Public entry point
 * ────────────────────────────────────────────────────────────────────────── */
int freesasa_shrake_rupley(double *sasa,
                           const coord_t *xyz,
                           const double *r,
                           const freesasa_parameters *param)
{
    int i, n_atoms, n_threads, resolution, return_value;
    double probe_radius;
    sr_data sr;

    assert(sasa); assert(xyz); assert(r);

    if (param == NULL) param = &freesasa_default_parameters;

    n_atoms      = freesasa_coord_n(xyz);
    n_threads    = param->n_threads;
    resolution   = param->shrake_rupley_n_points;
    probe_radius = param->probe_radius;
    return_value = FREESASA_SUCCESS;

    if (resolution <= 0)
        return fail_msg("%f test points invalid resolution in S&R, must be > 0\n", resolution);
    if (n_atoms == 0)
        return freesasa_warn("in %s(): empty coordinates", __func__);
    if (n_threads > n_atoms) {
        n_threads = n_atoms;
        freesasa_warn("no sense in having more threads than atoms, only using %d threads", n_threads);
    }

    if (init_sr(&sr, sasa, xyz, r, probe_radius, resolution, n_threads))
        return FREESASA_FAIL;

#if USE_OPENMP
    if (n_threads > 1) {
        #pragma omp parallel for schedule(dynamic, 64) num_threads(n_threads) \
            default(none) shared(sr, sasa, n_atoms)
        for (i = 0; i < n_atoms; ++i) {
            int tid = omp_get_thread_num();
            sasa[i] = sr_atom_area(i, &sr, tid);
        }
    } else {
        for (i = 0; i < n_atoms; ++i)
            sasa[i] = sr_atom_area(i, &sr, 0);
    }
#else
    if (n_threads > 1)
        freesasa_warn("in %s(): program compiled for single-threaded use, "
                      "but multiple threads were requested, will "
                      "proceed in single-threaded mode\n", __func__);
    for (i = 0; i < n_atoms; ++i)
        sasa[i] = sr_atom_area(i, &sr, 0);
#endif

    release_sr(&sr);
    return return_value;
}

/* ──────────────────────────────────────────────────────────────────────────
 * SIMD helpers
 *
 * check_neighbors_hidden():
 *   Returns 1 if test point (tx,ty,tz) is OUTSIDE ALL neighbor spheres
 *   (i.e. it is surface-exposed), 0 otherwise.
 *
 * We process 16 neighbors per iteration (AVX-512 float) or 8 (AVX2).
 * The loop is unrolled by the compiler with -O3 -march=native.
 * ────────────────────────────────────────────────────────────────────────── */

#if defined(FREESASA_USE_AVX512)

static inline int
check_neighbors_hidden(float tx, float ty, float tz,
                       const float * restrict cx,
                       const float * restrict cy,
                       const float * restrict cz,
                       const float * restrict cr2,
                       int n)
{
    /* Process 16 neighbors per iteration with AVX-512 */
    __m512 vtx = _mm512_set1_ps(tx);
    __m512 vty = _mm512_set1_ps(ty);
    __m512 vtz = _mm512_set1_ps(tz);
    int k = 0;

    for (; k + 16 <= n; k += 16) {
        __m512 dx = _mm512_sub_ps(_mm512_load_ps(cx + k), vtx);
        __m512 dy = _mm512_sub_ps(_mm512_load_ps(cy + k), vty);
        __m512 dz = _mm512_sub_ps(_mm512_load_ps(cz + k), vtz);
        __m512 d2 = _mm512_fmadd_ps(dx, dx,
                    _mm512_fmadd_ps(dy, dy,
                    _mm512_mul_ps(dz, dz)));
        __m512 r2v = _mm512_load_ps(cr2 + k);
        /* mask: lanes where d2 <= r2 (point is inside this neighbor) */
        __mmask16 inside = _mm512_cmp_ps_mask(d2, r2v, _CMP_LE_OQ);
        if (inside) return 0; /* hidden by at least one neighbor */
    }
    /* Scalar tail */
    for (; k < n; ++k) {
        float dx = cx[k] - tx, dy = cy[k] - ty, dz = cz[k] - tz;
        if (dx*dx + dy*dy + dz*dz <= cr2[k]) return 0;
    }
    return 1; /* surface-exposed */
}

#elif defined(FREESASA_USE_AVX2)

static inline int
check_neighbors_hidden(float tx, float ty, float tz,
                       const float * restrict cx,
                       const float * restrict cy,
                       const float * restrict cz,
                       const float * restrict cr2,
                       int n)
{
    __m256 vtx = _mm256_set1_ps(tx);
    __m256 vty = _mm256_set1_ps(ty);
    __m256 vtz = _mm256_set1_ps(tz);
    int k = 0;

    for (; k + 8 <= n; k += 8) {
        __m256 dx = _mm256_sub_ps(_mm256_load_ps(cx + k), vtx);
        __m256 dy = _mm256_sub_ps(_mm256_load_ps(cy + k), vty);
        __m256 dz = _mm256_sub_ps(_mm256_load_ps(cz + k), vtz);
        __m256 d2 = _mm256_fmadd_ps(dx, dx,
                    _mm256_fmadd_ps(dy, dy,
                    _mm256_mul_ps(dz, dz)));
        /* _CMP_LE_OQ = 0x12 */
        __m256 cmp = _mm256_cmp_ps(d2, _mm256_load_ps(cr2 + k), 0x12);
        if (_mm256_movemask_ps(cmp)) return 0;
    }
    for (; k < n; ++k) {
        float dx = cx[k] - tx, dy = cy[k] - ty, dz = cz[k] - tz;
        if (dx*dx + dy*dy + dz*dz <= cr2[k]) return 0;
    }
    return 1;
}

#else

/* Portable fallback: #pragma omp simd for auto-vectorization */
static inline int
check_neighbors_hidden(float tx, float ty, float tz,
                       const float * restrict cx,
                       const float * restrict cy,
                       const float * restrict cz,
                       const float * restrict cr2,
                       int n)
{
    int exposed = 1;
    #pragma omp simd reduction(&:exposed)
    for (int k = 0; k < n; ++k) {
        float dx = cx[k] - tx, dy = cy[k] - ty, dz = cz[k] - tz;
        if (dx*dx + dy*dy + dz*dz <= cr2[k]) exposed = 0;
    }
    return exposed;
}

#endif /* SIMD selection */

/* ──────────────────────────────────────────────────────────────────────────
 * Per-atom SASA contribution
 *
 * Uses per-atom SoA float32 neighbor cache for SIMD inner loop.
 * Outer j-loop over test points uses the original "sticky neighbor" trick
 * for cache locality but now the inner k-scan is fully vectorized.
 * ────────────────────────────────────────────────────────────────────────── */
static double
sr_atom_area(int i, const sr_data *sr, int thread_index)
{
    const int n_points   = sr->n_points;
    int *spcount         = sr->spcount[thread_index];
    const nb_cache_t *c  = &sr->nb_cache[i];
    const int nni        = c->n;
    const float * restrict cx  = c->x;
    const float * restrict cy  = c->y;
    const float * restrict cz  = c->z;
    const float * restrict cr2 = c->r2;
    const double ri      = sr->r[i];
    const double *restrict v = freesasa_coord_all(sr->xyz);
    const double *restrict vi = v + 3 * i;

    /* Scale and translate test points for atom i */
    coord_t *restrict tp_coord = sr->tp_local[thread_index];
    freesasa_coord_copy(tp_coord, sr->srp);
    freesasa_coord_scale(tp_coord, ri);
    freesasa_coord_translate(tp_coord, vi);
    const double *restrict tp = freesasa_coord_all(tp_coord);

    memset(spcount, 0, n_points * sizeof(int));

    if (nni == 0) {
        /* No neighbors — all test points are exposed */
        return 4.0 * M_PI * ri * ri;
    }

    /* ── Main loop over test points ──────────────────────────────────────
     * "Sticky neighbor" trick: remember the last occluding neighbor index
     * and check it first. On a miss, fall back to full SIMD scan.
     * This preserves spatial locality (consecutive test points on the spiral
     * tend to be occluded by the same neighbor).
     * ──────────────────────────────────────────────────────────────────── */
    int current_nb = 0;  /* sticky neighbor index into SoA cache */
    int n_surface = 0;

    for (int j = 0; j < n_points; ++j) {
        float tx = (float)tp[j * 3];
        float ty = (float)tp[j * 3 + 1];
        float tz = (float)tp[j * 3 + 2];

        /* Fast path: check sticky neighbor first */
        {
            float dx = cx[current_nb] - tx;
            float dy = cy[current_nb] - ty;
            float dz = cz[current_nb] - tz;
            if (dx*dx + dy*dy + dz*dz <= cr2[current_nb]) {
                /* Hidden by sticky neighbor — move to next point */
                continue;
            }
        }

        /* Slow path: full SIMD scan over all nni neighbors */
        if (check_neighbors_hidden(tx, ty, tz, cx, cy, cz, cr2, nni)) {
            spcount[j] = 1;  /* surface exposed */
        } else {
            /* Find which neighbor occluded (update sticky) */
            for (int k = 0; k < nni; ++k) {
                float dx = cx[k] - tx, dy = cy[k] - ty, dz = cz[k] - tz;
                if (dx*dx + dy*dy + dz*dz <= cr2[k]) {
                    current_nb = k;
                    break;
                }
            }
        }
    }

    for (int k = 0; k < n_points; ++k)
        if (spcount[k]) ++n_surface;

    return (4.0 * M_PI * ri * ri * n_surface) / n_points;
}
