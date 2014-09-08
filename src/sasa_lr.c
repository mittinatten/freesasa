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

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#if HAVE_CONFIG_H
# include <config.h>
#endif
#if HAVE_PTHREAD_H
# include <pthread.h>
#endif

#include "freesasa.h"
#include "sasa.h"
#include "srp.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

extern const char *freesasa_name;
extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

//calculation parameters (results stored in *sasa)
typedef struct {
    int n_atoms;
    double *radii; //including probe
    const freesasa_coord_t *xyz;
    int **nb; // neighbors
    int *nn; // number of neighbors
    double **nb_xyd; // neighbours, xy-distance
    double **nb_xd; // neighbours, x-distance
    double **nb_yd; // neighbours, y-distance
    double delta; // slice width
    double min_z; // bounds of the molecule
    double max_z;
    double *sasa; // results
} sasa_lr_t;

typedef struct {
    int n_slice; //number of atoms in slice
    double z; //the mid-point of the slice
    char *in_slice;
    int *idx; //index in slice to global numbering
    int *xdi; //global numbering to index in slice
    double *DR; //corrective multiplicative factor (D in L&R paper)
    double *r; //radius in slice;
    double *exposed_arc; //exposed arc length (in Å) for each atom
} sasa_lr_slice_t;

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr_t*);
static void *sasa_lr_thread(void *arg);
#endif

/** Init slice and sum up arc-length */
static void sasa_add_slice_area(double z, sasa_lr_t*);

/** Find which arcs are exposed in a slice */
static void sasa_exposed_arcs(const sasa_lr_slice_t *,
                              const sasa_lr_t *);

/** a and b are a set of alpha and betas (in the notation of the
    manual). This function finds the union of those intervals on the
    circle, and returns 2*PI minus the length of the joined
    interval(s) (i.e. the exposed arc length). Does not necessarily
    leave a and b in a consistent state. */
static double sasa_sum_angles(int n_buried, double * restrict a,
                              double * restrict b);

/** Calculate adjacency lists, given coordinates and radii. Also
    stores some distances to be used later. */
static void sasa_get_contacts(sasa_lr_t *lr);

/** Initialize object to be used for L&R calculation */
static sasa_lr_t* freesasa_init_lr(double *sasa,
                                   const freesasa_coord_t *xyz,
                                   const double *atom_radii,
                                   double probe_radius,
                                   double delta) {
    sasa_lr_t* lr = (sasa_lr_t*) malloc(sizeof(sasa_lr_t));
    const int n_atoms = freesasa_coord_n(xyz);
    double max_z=-1e50, min_z=1e50;
    double max_r = 0;
    double *radii = (double*) malloc(sizeof(double)*n_atoms);
    const double *v = freesasa_coord_all(xyz);
    //find bounds of protein along z-axis and init radii
    for (size_t i = 0; i < n_atoms; ++i) {
        radii[i] = atom_radii[i] + probe_radius;
        double z = v[3*i+2], r = radii[i];
        max_z = z > max_z ? z : max_z;
        min_z = z < min_z ? z : min_z;
        sasa[i] = 0.;
        max_r = r > max_r ? r : max_r;
    }
    min_z -= max_r;
    max_z += max_r;
    min_z += 0.5*delta;

    //copy parameters
    lr->n_atoms = n_atoms; lr->radii = radii; lr->xyz = xyz;
    lr->delta = delta;
    lr->min_z = min_z; lr->max_z = max_z;
    lr->sasa = sasa;

    //these will be malloc'd by get_contacts
    lr->nn = NULL;
    lr->nb = NULL;
    lr->nb_xyd = lr->nb_xd = lr->nb_yd = NULL;

    return lr;
}

static void freesasa_free_lr(sasa_lr_t *lr)
{
    for (int i = 0; i < lr->n_atoms; ++i) {
        free(lr->nb[i]);
        free(lr->nb_xyd[i]);
        free(lr->nb_xd[i]);
        free(lr->nb_yd[i]);
    }
    free(lr->radii);
    free(lr->nn);
    free(lr->nb);
    free(lr->nb_xyd);
    free(lr->nb_xd);
    free(lr->nb_yd);
    free(lr);
}

int freesasa_lee_richards(double *sasa,
                          const freesasa_coord_t *xyz,
                          const double *atom_radii,
                          double probe_radius,
                          double delta,
                          int n_threads)
{
    int return_value = FREESASA_SUCCESS;
    if (freesasa_coord_n(xyz) == 0) {
        return freesasa_warn("Attempting Lee & Richards calculation "
                             "on empty coordinates");
    }

    // determine slice range and init radii and sasa arrays
    sasa_lr_t *lr = freesasa_init_lr(sasa, xyz, atom_radii,
                                     probe_radius, delta);

    // determine which atoms are neighbours
    sasa_get_contacts(lr);

    if (n_threads > 1) {
#if HAVE_LIBPTHREAD
        sasa_lr_do_threads(n_threads, lr);
#else
        return_value = freesasa_warn("program compiled for single-threaded use, "
                                     "but multiple threads were requested. Will "
                                     "proceed in single-threaded mode.\n");
        n_threads = 1;
#endif
    }
    if (n_threads == 1) {
        // loop over slices
        for (double z = lr->min_z; z < lr->max_z; z += lr->delta) {
            sasa_add_slice_area(z,lr);
        }
    }
    freesasa_free_lr(lr);

    return return_value;
}

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr_t *lr)
{
    double *t_sasa[n_threads];
    pthread_t thread[n_threads];
    sasa_lr_t lrt[n_threads];
    const double max_z = lr->max_z, min_z = lr->min_z, delta = lr->delta;
    int n_slices = (int)ceil((max_z-min_z)/delta);
    int n_perthread = n_slices/n_threads;
    for (int t = 0; t < n_threads; ++t) {
        t_sasa[t] = (double*)malloc(sizeof(double)*lr->n_atoms);
        for (int i = 0; i < lr->n_atoms; ++i) t_sasa[t][i] = 0;
        lrt[t] = *lr;
        lrt[t].sasa = t_sasa[t];
        lrt[t].min_z = min_z + t*n_perthread*delta;
        if (t < n_threads - 1) {
            lrt[t].max_z = lrt[t].min_z + n_perthread*delta;
        } else {
            lrt[t].max_z = max_z;
        }
        int res = pthread_create(&thread[t], NULL, sasa_lr_thread, (void *) &lrt[t]);
        if (res) {
            perror(freesasa_name);
            exit(EXIT_FAILURE);
        }
    }
    for (int t = 0; t < n_threads; ++t) {
        void *thread_result;
        int res = pthread_join(thread[t],&thread_result);
        if (res) {
            perror(freesasa_name);
            exit(EXIT_FAILURE);
        }
    }

    for (int t = 0; t < n_threads; ++t) {
        for (int i = 0; i < lr->n_atoms; ++i) {
            lr->sasa[i] += t_sasa[t][i];
        }
        free(t_sasa[t]);
    }
}

static void *sasa_lr_thread(void *arg)
{
    sasa_lr_t *lr = ((sasa_lr_t*) arg);
    for (double z = lr->min_z; z < lr->max_z; z += lr->delta) {
        sasa_add_slice_area(z, lr);
    }
    pthread_exit(NULL);
}
#endif

// find which atoms are in slice
static sasa_lr_slice_t* sasa_init_slice(double z, const sasa_lr_t *lr)
{
    sasa_lr_slice_t *slice = (sasa_lr_slice_t*) malloc(sizeof(sasa_lr_slice_t));
    int n_atoms = lr->n_atoms;
    int n_slice = slice->n_slice = 0;
    slice->xdi = (int*) malloc(n_atoms*sizeof(int));
    slice->in_slice = (char*) malloc(n_atoms*sizeof(int));
    double delta = lr->delta;
    const double *v = freesasa_coord_all(lr->xyz);
    slice->idx = NULL;
    slice->r = slice->DR = NULL;
    slice->z = z;

    memset(slice->in_slice,0,sizeof(int)*n_atoms);

    // locate atoms in each slice and do some initialization
    for (size_t i = 0; i < n_atoms; ++i) {
        double ri = lr->radii[i];
        double d = fabs(v[3*i+2]-z);
        double r;
        if (d < ri) {
            slice->idx = (int*) realloc(slice->idx,(n_slice+1)*sizeof(int));
            slice->r = (double*) realloc(slice->r,(n_slice+1)*sizeof(double));
            slice->DR = (double*) realloc(slice->DR,(n_slice+1)*sizeof(double));

            slice->idx[n_slice] = i;
            slice->xdi[i] = n_slice;
            slice->in_slice[i] = 1;
            slice->r[n_slice] = r = sqrt(ri*ri-d*d);
            slice->DR[n_slice] = ri/r*(delta/2. +
                                       (delta/2. < ri-d ? delta/2. : ri-d));
            ++n_slice;
        }  else {
            slice->xdi[i] = -1;
        } 
    }
    slice->n_slice = n_slice;
    slice->exposed_arc = (double* )malloc(n_slice*(sizeof(double)));
    return slice;
}

static void sasa_free_slice(sasa_lr_slice_t* slice)
{
    free(slice->idx);
    free(slice->xdi);
    free(slice->in_slice);
    free(slice->r);
    free(slice->DR);
    free(slice->exposed_arc);
    free(slice);
}

static void sasa_add_slice_area(double z, sasa_lr_t *lr)
{
    sasa_lr_slice_t *slice = sasa_init_slice(z,lr);

    //find exposed arcs
    sasa_exposed_arcs(slice, lr);

    // calculate contribution to each atom's SASA from the present slice
    for (int i = 0; i < slice->n_slice; ++i) {
        lr->sasa[slice->idx[i]] += slice->exposed_arc[i];
    }
    sasa_free_slice(slice);
}

static void sasa_exposed_arcs(const sasa_lr_slice_t *slice,
                              const sasa_lr_t *lr)
{
    const int n_slice = slice->n_slice;
    const int *nn = lr->nn;
    int * const *nb = lr->nb;
    double * const *nb_xyd = lr->nb_xyd;
    const double *r = slice->r;
    const double *v = freesasa_coord_all(lr->xyz);

    char is_completely_buried[n_slice]; // keep track of completely buried circles
    memset(is_completely_buried,0,sizeof is_completely_buried);
    
    //loop over atoms in slice
    //lower-case i,j is atoms in the slice, upper-case I,J are their
    //corresponding global indexes
    for (int i = 0; i < n_slice; ++i) {
        slice->exposed_arc[i] = 0;
        if (is_completely_buried[i]) {
            continue;
        }
        int I = slice->idx[i];
        double ri = slice->r[i], a[nn[I]], b[nn[I]];
        int n_buried = 0;
        // loop over neighbors
        for (int ni = 0; ni < nn[I]; ++ni) {
            int J = nb[I][ni];
            if (slice->in_slice[J] == 0) continue;
            int j = slice->xdi[J];
            double rj = r[j];
            double d = nb_xyd[I][ni];
            // reasons to skip calculation
            if (d >= ri + rj) continue;     // atoms aren't in contact
            if (d + ri < rj) { // circle i is completely inside j
                is_completely_buried[i] = 1;
                break;
            }
            if (d + rj < ri) { // circle j is completely inside i
                is_completely_buried[j] = 1;
                continue;
            }
            // half the arclength occluded from circle i due to verlap with circle j
            double alpha = acos ((ri*ri + d*d - rj*rj)/(2.0*ri*d));
            // the polar coordinates angle of the vector connecting i and j
            double beta = atan2 (lr->nb_yd[I][ni],lr->nb_xd[I][ni]);
            a[n_buried] = alpha;
            b[n_buried] = beta;

            ++n_buried;
        }
        if (is_completely_buried[i] == 0) {
            slice->exposed_arc[i] = 
                ri*slice->DR[i]*sasa_sum_angles(n_buried,a,b);
        }
#ifdef DEBUG
        if (is_completely_buried[i] == 0) {
            //exposed_arc[i] = 0;
            const double *v = freesasa_coord_all(lr->xyz);
            double xi = v[3*I], yi = v[3*I+1];
            for (double c = 0; c < 2*PI; c += PI/45.0) {
                int is_exp = 1;
                for (int i = 0; i < n_buried; ++i) {
                    if ((c > b[i]-a[i] && c < b[i]+a[i]) ||
                        (c - 2*PI > b[i]-a[i] && c - 2*PI < b[i]+a[i]) ||
                        (c + 2*PI > b[i]-a[i] && c + 2*PI < b[i]+a[i])) {
                        is_exp = 0; break;
                    }
                }
                // print the arcs used in calculation
                if (is_exp) printf("%6.2f %6.2f %6.2f %7.5f\n",
                                   xi+ri*cos(c),yi+ri*sin(c),slice->z,c);
            }
            printf("\n");
        }
#endif
    }
}

static double sasa_sum_angles(int n_buried, double *restrict a, double *restrict b)
{
    /* Innermost function in L&R, could potentially be sped up, but
       probably requires rethinking, algorithmically. Perhaps
       recursion could be rolled out somehow. */
    char excluded[n_buried], n_exc = 0, n_overlap = 0;
    memset(excluded,0,n_buried);

    for (int i = 0; i < n_buried; ++i) {
        if (excluded[i]) continue;
        for (int j = 0; j < n_buried; ++j) {
            if (excluded[j]) continue;
            if (i == j) continue;

            //check for overlap
            double bi = b[i], ai = a[i]; //will be updating throughout the loop
            double bj = b[j], aj = a[j];
            double d;
            for (;;) {
                d = bj - bi;
                if (d > PI) bj -= 2*PI;
                else if (d < -PI) bj += 2*PI;
                else break;
            }
            if (fabs(d) > ai+aj) continue;
            ++n_overlap;

            //calculate new joint interval
            double inf_i = bi-ai, inf_j = bj-aj;
            double sup_i = bi+ai, sup_j = bj+aj;
            double inf = inf_i < inf_j ? inf_i : inf_j;
            double sup = sup_i > sup_j ? sup_i : sup_j;
            b[i] = (inf + sup)/2.0;
            a[i] = (sup - inf)/2.0;
            if (a[i] > PI) return 0;
            if (b[i] > PI) b[i] -= 2*PI;
            if (b[i] < -PI) b[i] += 2*PI;

            a[j] = 0; // the j:th interval should be ignored
            excluded[j] = 1;
            if (++n_exc == n_buried-1) break;
        }
        if (n_exc == n_buried-1) break; // means everything's been counted
    }

    // recursion until no overlapping intervals
    if (n_overlap) {
        double b2[n_buried], a2[n_buried];
        int n = 0;
        for (int i = 0; i < n_buried; ++i) {
            if (excluded[i] == 0) {
                b2[n] = b[i];
                a2[n] = a[i];
                ++n;
            }
        }
        return sasa_sum_angles(n,a2,b2);
    }
    // else return angle
    double buried_angle = 0;
    for (int i = 0; i < n_buried; ++i) {
        buried_angle += 2.0*a[i];
    }
    return 2*PI - buried_angle;
}
static void sasa_get_contacts(sasa_lr_t *lr)
{
    /* For low resolution L&R this function is the bottleneck in
       speed. Will also depend on number of atoms. */
    size_t n_atoms = lr->n_atoms;
    const double *v = freesasa_coord_all(lr->xyz);
    const double *radii = lr->radii;

    // init adjacency lists
    int *nn = (int*)malloc(sizeof(int)*n_atoms);
    int **nb = (int**)malloc(sizeof(int*)*n_atoms);
    double **nb_xyd = (double**) malloc(sizeof(double*)*n_atoms);
    double **nb_xd = (double**) malloc(sizeof(double*)*n_atoms);
    double **nb_yd = (double**) malloc(sizeof(double*)*n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        nn[i] = 0; nb[i] = NULL;
        nb_xyd[i] = NULL; nb_xd[i] = NULL; nb_yd[i] = NULL;
    }

    //calculate lists
    for (int i = 0; i < n_atoms; ++i) {
        double ri = radii[i];
        double xi = v[i*3], yi = v[i*3+1], zi = v[i*3+2];
        for (int j = i+1; j < n_atoms; ++j) {
            double rj = radii[j];
            double cut2 = (ri+rj)*(ri+rj);

            /* most pairs of atoms will be far away from each other on
               at least one axis, the following improves speed
               significantly for large proteins */
            double xj = v[j*3], yj = v[j*3+1], zj = v[j*3+2];
            if ((xj-xi)*(xj-xi) > cut2 ||
                (yj-yi)*(yj-yi) > cut2 ||
                (zj-zi)*(zj-zi) > cut2) {
                continue;
            }
            double dx = xj-xi, dy = yj-yi, dz = zj-zi;
            if (dx*dx + dy*dy + dz*dz < cut2) {
                ++nn[i]; ++nn[j];
                //record neighbors
                nb[i] = (int*) realloc(nb[i],sizeof(int)*nn[i]);
                nb[j] = (int*) realloc(nb[j],sizeof(int)*nn[j]);
                nb[i][nn[i]-1] = j;
                nb[j][nn[j]-1] = i;

                //record neighbor distances
                nb_xyd[i] = (double*) realloc(nb_xyd[i],sizeof(double)*nn[i]);
                nb_xyd[j] = (double*) realloc(nb_xyd[j],sizeof(double)*nn[j]);
                nb_xd[i] = (double*) realloc(nb_xd[i],sizeof(double)*nn[i]);
                nb_xd[j] = (double*) realloc(nb_xd[j],sizeof(double)*nn[j]);
                nb_yd[i] = (double*) realloc(nb_yd[i],sizeof(double)*nn[i]);
                nb_yd[j] = (double*) realloc(nb_yd[j],sizeof(double)*nn[j]);

                double d = sqrt(dx*dx+dy*dy);
                nb_xyd[i][nn[i]-1] = d;
                nb_xyd[j][nn[j]-1] = d;

                nb_xd[i][nn[i]-1] = dx;
                nb_xd[j][nn[j]-1] = -dx;
                nb_yd[i][nn[i]-1] = dy;
                nb_yd[j][nn[j]-1] = -dy;
            }
        }
    }

    //copy results
    lr->nn = nn;
    lr->nb = nb;
    lr->nb_xyd = nb_xyd;
    lr->nb_xd = nb_xd;
    lr->nb_yd = nb_yd;
}

