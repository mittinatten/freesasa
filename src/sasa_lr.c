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
#include "verlet.h"

extern const char *freesasa_name;
extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

const double TWOPI = 2*M_PI;

//calculation parameters and data (results stored in *sasa)
typedef struct {
    int n_atoms;
    double *radii; //including probe
    const freesasa_coord *xyz;
    freesasa_verlet *adj;
    double delta; // slice width
    double min_z; // bounds of the molecule
    double max_z;
    double *sasa; // results
} lr_data;

typedef struct {
    int first_atom;
    int last_atom;
    lr_data *lr;
} lr_thread_interval;

#if HAVE_LIBPTHREAD
static void lr_do_threads(int n_threads, lr_data*);
static void *lr_thread(void *arg);
#endif

/** Returns the are of atom i */
static double atom_area(lr_data *lr,int i);

/** Sum of exposed arcs based on buried arc intervals arc, assumes no
    intervals cross zero */
static double sum_arcs(int n_buried, double *restrict arc);

/** Initialize object to be used for L&R calculation */
static lr_data* init_lr(double *sasa,
                        const freesasa_coord *xyz,
                        const double *atom_radii,
                        double probe_radius,
                        double delta) {
    const int n_atoms = freesasa_coord_n(xyz);
    double max_z=-1e50, min_z=1e50;
    double max_r = 0;
    const double *v = freesasa_coord_all(xyz);
    lr_data* lr = malloc(sizeof(lr_data));
    double *radii = malloc(sizeof(double)*n_atoms);
    assert(lr);
    assert(radii);
    //find bounds of protein along z-axis and init radii
    for (int i = 0; i < n_atoms; ++i) {
        double z, r;
        radii[i] = atom_radii[i] + probe_radius;
        z = v[3*i+2];
        r = radii[i];
        max_z = fmax(z,max_z);
        min_z = fmin(z,min_z);
        max_r = fmax(r, max_r);
        sasa[i] = 0.;
    }
    min_z -= max_r;
    max_z += max_r;
    min_z += 0.5*delta;

    //copy parameters
    lr->n_atoms = n_atoms; lr->radii = radii; lr->xyz = xyz;
    lr->delta = delta;
    lr->min_z = min_z; lr->max_z = max_z;
    lr->sasa = sasa;
    lr->adj = NULL;

    return lr;
}

static void free_lr(lr_data *lr)
{
    free(lr->radii);
    freesasa_verlet_free(lr->adj);
    free(lr);
}

int freesasa_lee_richards(double *sasa,
                          const freesasa_coord *xyz,
                          const double *atom_radii,
                          double probe_radius,
                          double delta,
                          int n_threads)
{
    assert(sasa);
    assert(xyz);
    assert(atom_radii);

    int return_value = FREESASA_SUCCESS;
    lr_data *lr;

    if (freesasa_coord_n(xyz) == 0) {
        return freesasa_warn("%s: empty coordinates",__func__);
    }

    // determine slice range and init radii and sasa arrays
    lr = init_lr(sasa, xyz, atom_radii, probe_radius, delta);

    // determine which atoms are neighbours
    lr->adj = freesasa_verlet_new(xyz,lr->radii);

    if (n_threads > 1) {
#if HAVE_LIBPTHREAD
        lr_do_threads(n_threads, lr);
#else
        return_value = freesasa_warn("%s: program compiled for single-threaded use, "
                                     "but multiple threads were requested. Will "
                                     "proceed in single-threaded mode.\n",
                                     __func__);
        n_threads = 1;
#endif /* pthread */
    }
    if (n_threads == 1) {
        for (int i = 0; i < lr->n_atoms; ++i) {
            lr->sasa[i] = atom_area(lr,i);
        }        
    }
    free_lr(lr);

    return return_value;
}

#if HAVE_LIBPTHREAD
static void lr_do_threads(int n_threads, lr_data *lr)
{
    pthread_t thread[n_threads];
    lr_thread_interval t_data[n_threads];
    int n_atoms = lr->n_atoms, n_perthread = n_atoms/n_threads, res;
    void *thread_result;
 
    for (int t = 0; t < n_threads; ++t) {
        t_data[t].first_atom = t*n_perthread;
        if (t == n_threads-1) {
            t_data[t].last_atom = n_atoms - 1;
        } else {
            t_data[t].last_atom = (t+1)*n_perthread - 1;
        }
        t_data[t].lr = lr;
        res = pthread_create(&thread[t], NULL, lr_thread,
                                 (void *) &t_data[t]);
        if (res) {
            perror(freesasa_name);
            exit(EXIT_FAILURE);
        }
    }
    for (int t = 0; t < n_threads; ++t) {
        res = pthread_join(thread[t],&thread_result);
        if (res) {
            perror(freesasa_name);
            exit(EXIT_FAILURE);
        }
    }
}

static void *lr_thread(void *arg)
{
    lr_thread_interval *ti = ((lr_thread_interval*) arg);
    for (int i = ti->first_atom; i <= ti->last_atom; ++i) {
        /* the different threads write to different parts of the
           array, so locking shouldn't be necessary */
        ti->lr->sasa[i] = atom_area(ti->lr, i);
    }
    pthread_exit(NULL);
}
#endif /* pthread */

static double atom_area(lr_data *lr,int i)
{
    // This function is large because a large number of pre-calculated
    // arrays need to be accessed efficiently. Partially dereferenced
    // here to make access more efficient.
    const int nni = lr->adj->nn[i];
    const double * restrict const v = freesasa_coord_all(lr->xyz);
    const double * restrict const r = lr->radii;
    const int * restrict const nbi = lr->adj->nb[i];
    const double * restrict const xydi = lr->adj->nb_xyd[i];
    const double * restrict const xdi = lr->adj->nb_xd[i];
    const double * restrict const ydi = lr->adj->nb_yd[i];
    const double zi = v[3*i+2], delta = lr->delta, ri = r[i], d_half = delta/2.;
    double arc[nni*4], z_nb[nni], r_nb[nni];
    double z_slice, z0, sasa = 0;
    const int bottom = ((zi-ri)-lr->min_z)/delta + 1;

    z0 = lr->min_z+bottom*delta;
    
    for (int j = 0; j < nni; ++j) {
        z_nb[j] = v[3*nbi[j]+2];
        r_nb[j] = r[nbi[j]];
    }
    
    for (z_slice = z0; z_slice < zi+ri; z_slice += delta) {
        const double di = fabs(zi - z_slice);
        const double ri_slice2 = ri*ri-di*di;
        if (ri_slice2 < 0 ) continue; // handle round-off errors
        const double ri_slice = sqrt(ri_slice2);
        const double DR = ri/ri_slice*(d_half + fmin(d_half,ri-di));
        int n_arcs = 0, is_buried = 0;
        for (int j = 0; j < nni; ++j) {
            const double zj = z_nb[j];
            const double dj = fabs(zj - z_slice);
            const double rj = r_nb[j];
            if (dj < rj) {
                const double rj_slice2 = rj*rj-dj*dj;
                const double rj_slice = sqrt(rj_slice2);
                const double dij = xydi[j];
                double alpha, beta, inf, sup;
                int narc2;
                if (dij >= ri_slice + rj_slice) { // atoms aren't in contact
                    continue;
                }
                if (dij + ri_slice < rj_slice) { // circle i is completely inside j
                    is_buried = 1;
                    break;
                }
                if (dij + rj_slice < ri_slice) { // circle j is completely inside i
                    continue;
                }
                // arc of circle i intersected by circle j
                alpha = acos ((ri_slice2 + dij*dij - rj_slice2)/(2.0*ri_slice*dij));
                // position of mid-point of intersection along circle i
                beta = atan2 (ydi[j],xdi[j]) + M_PI;
                inf = beta - alpha;
                sup = beta + alpha;
                if (inf < 0) inf += TWOPI;
                if (sup > 2*M_PI) sup -= TWOPI;
                narc2 = 2*n_arcs;
                // store the arc, if arc passes 2*PI split into two
                if (sup < inf) {
                    //store arcs as contiguous pairs of angles
                    arc[narc2]   = 0;
                    arc[narc2+1] = sup;
                    //second arc
                    arc[narc2+2] = inf;
                    arc[narc2+3] = TWOPI;
                    n_arcs += 2;
                } else { 
                    arc[narc2]   = inf;
                    arc[narc2+1] = sup;
                    ++n_arcs;
                }
            }
        }
        if (is_buried == 0) {
            sasa += ri_slice*DR*sum_arcs(n_arcs,arc);
        }
#ifdef DEBUG
        if (completely_buried == 0) {
            //exposed_arc[i] = 0;
            double xi = v[3*i], yi = v[3*i+1];
            for (double c = 0; c < 2*M_PI; c += M_PI/30.0) {
                int is_exp = 1;
                for (int j = 0; j < n_buried; ++j) {
                    double inf = arc[2*j], sup = arc[2*j+1];
                    if (c >= inf && c <= sup) {
                        is_exp = 0; break;
                    }
                }
                double d = c-M_PI;
                // print the arcs used in calculation
                if (is_exp) printf("%6.2f %6.2f %6.2f %7.5f\n",
                                   xi+ri_slice*cos(d),yi+ri_slice*sin(d),z_slice,d);
            }
            printf("\n");
        }
#endif /* Debug */

    }
    return sasa;
}

//insertion sort (faster than qsort for these short lists)
inline static void sort_arcs(int n, double *restrict arc) 
{
    double tmp[2];
    double *end = arc+2*n, *arcj, *arci;
    for (arci = arc+2; arci < end; arci += 2) {
        memcpy(tmp,arci,2*sizeof(double));
        arcj = arci;
        while (arcj > arc && *(arcj-2) > tmp[0]) {
            memcpy(arcj,arcj-2,2*sizeof(double));
            arcj -= 2;
        }
        memcpy(arcj,tmp,2*sizeof(double));
    }
}

// sort arcs by start-point, loop through them to sum parts of circle
// not covered by any of the arcs
inline static double sum_arcs(int n, double *arc)
{
    if (n == 0) return TWOPI;
    double sum, sup, tmp;
    sort_arcs(n,arc);
    sum = arc[0];
    sup = arc[1];
    for (int i2 = 2; i2 < 2*n; i2 += 2) {
        if (sup < arc[i2]) sum += arc[i2] - sup;
        tmp = arc[i2+1];
        if (tmp > sup) sup = tmp;
    } 
    return sum + TWOPI - sup;
}
