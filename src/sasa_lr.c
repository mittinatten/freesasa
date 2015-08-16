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
#include "adjacency.h"

extern const char *freesasa_name;
extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

//calculation parameters (results stored in *sasa)
typedef struct {
    int n_atoms;
    double *radii; //including probe
    const freesasa_coord *xyz;
    freesasa_adjacency *adj;
    double delta; // slice width
    double min_z; // bounds of the molecule
    double max_z;
    double *sasa; // results
} sasa_lr;

typedef struct {
    int first_atom;
    int last_atom;
    sasa_lr *lr;
} sasa_lr_thread_interval;

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr*);
static void *sasa_lr_thread(void *arg);
#endif

/** Returns the are of atom i */
static double sasa_atom_area(sasa_lr *lr,int i);

/** a and b are a set of alpha and betas (in the notation of the
    manual). This function finds the union of those intervals on the
    circle, and returns 2*PI minus the length of the joined
    interval(s) (i.e. the exposed arc length). Does not necessarily
    leave a and b in a consistent state. */
static double sasa_sum_angles(int n_buried, double * restrict a,
                              double * restrict b);

/** Same as the above but discretizes circle */
static double sasa_sum_angles_int(int n_buried, double *restrict a,
                                  double *restrict b);

/** Initialize object to be used for L&R calculation */
static sasa_lr* freesasa_init_lr(double *sasa,
                                 const freesasa_coord *xyz,
                                 const double *atom_radii,
                                 double probe_radius,
                                 double delta) {
    const int n_atoms = freesasa_coord_n(xyz);
    double max_z=-1e50, min_z=1e50;
    double max_r = 0;
    const double *v = freesasa_coord_all(xyz);
    sasa_lr* lr = malloc(sizeof(sasa_lr));
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
        sasa[i] = 0.;
        max_r = fmax(r, max_r);
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

static void freesasa_free_lr(sasa_lr *lr)
{
    free(lr->radii);
    freesasa_adjacency_free(lr->adj);
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
    sasa_lr *lr;

    if (freesasa_coord_n(xyz) == 0) {
        return freesasa_warn("%s: empty coordinates",__func__);
    }

    // determine slice range and init radii and sasa arrays
    lr = freesasa_init_lr(sasa, xyz, atom_radii,
                          probe_radius, delta);

    // determine which atoms are neighbours
    lr->adj = freesasa_adjacency_new(xyz,lr->radii);

    if (n_threads > 1) {
#if HAVE_LIBPTHREAD
        sasa_lr_do_threads(n_threads, lr);
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
            lr->sasa[i] = sasa_atom_area(lr,i); 
        }        
    }
    freesasa_free_lr(lr);

    return return_value;
}

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr *lr)
{
    pthread_t thread[n_threads];
    sasa_lr_thread_interval t_data[n_threads];
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
        res = pthread_create(&thread[t], NULL, sasa_lr_thread,
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

static void *sasa_lr_thread(void *arg)
{
    sasa_lr_thread_interval *ti = ((sasa_lr_thread_interval*) arg);
    for (int i = ti->first_atom; i <= ti->last_atom; ++i) {
        /* the different threads write to different parts of the
           array, so locking shouldn't be necessary */
        ti->lr->sasa[i] = sasa_atom_area(ti->lr, i);
    }
    pthread_exit(NULL);
}
#endif /* pthread */

static double sasa_atom_area(sasa_lr *lr,int i)
{
    int nni = lr->adj->nn[i];
    const double * restrict v = freesasa_coord_all(lr->xyz);
    const double * restrict r = lr->radii;
    const int * restrict nbi = lr->adj->nb[i];
    const double * restrict xydi = lr->adj->nb_xyd[i];
    const double * restrict xdi = lr->adj->nb_xd[i];
    const double * restrict ydi = lr->adj->nb_yd[i];
    double sasa = 0, zi = v[3*i+2], z_slice, z0,
        delta = lr->delta, ri = r[i], d_half = delta/2.;
    double a[nni], b[nni], z_nb[nni], r_nb[nni];
    int bottom = ((zi-ri)-lr->min_z)/delta + 1;

    z0 = lr->min_z+bottom*delta;
    
    for (int j = 0; j < nni; ++j) {
        z_nb[j] = v[3*nbi[j]+2];
        r_nb[j] = r[nbi[j]];
    }
    
    for (z_slice = z0; z_slice < zi+ri; z_slice += delta) {
        double di = fabs(zi - z_slice);
        double ri_slice2 = ri*ri-di*di;
        if (ri_slice2 < 0 ) continue; // handle round-off errors
        double ri_slice = sqrt(ri_slice2);
        double DR = ri/ri_slice*(d_half + fmin(d_half,ri-di));
        int n_buried = 0, completely_buried = 0;
        for (int j = 0; j < nni; ++j) {
            double zj = z_nb[j];
            double dj = fabs(zj - z_slice);
            double rj = r_nb[j];
            if (dj < rj) {
                double rj_slice = sqrt(rj*rj-dj*dj);
                double dij = xydi[j];
                if (dij >= ri_slice + rj_slice) { // atoms aren't in contact
                    continue;
                }
                if (dij + ri_slice < rj_slice) { // circle i is completely inside j
                    completely_buried = 1;
                    break;
                }
                if (dij + rj_slice < ri_slice) { // circle j is completely inside i
                    continue;
                }
                a[n_buried] = acos ((ri_slice*ri_slice + dij*dij
                                     - rj_slice*rj_slice)/(2.0*ri_slice*dij));
                b[n_buried] = atan2 (ydi[j],xdi[j]);
                ++n_buried;
            }
        }
        if (completely_buried == 0) {
            sasa += ri_slice*DR*sasa_sum_angles(n_buried,a,b);
        }
#ifdef DEBUG
        if (completely_buried == 0) {
            //exposed_arc[i] = 0;
            double xi = v[3*i], yi = v[3*i+1];
            for (double c = 0; c < 2*M_PI; c += M_PI/30.0) {
                int is_exp = 1;
                for (int j = 0; j < n_buried; ++j) {
                    if ((c > b[j]-a[j] && c < b[j]+a[j]) ||
                        (c - 2*M_PI > b[j]-a[j] && c - 2*M_PI < b[j]+a[j]) ||
                        (c + 2*M_PI > b[j]-a[j] && c + 2*M_PI < b[j]+a[j])) {
                        is_exp = 0; break;
                    }
                }
                // print the arcs used in calculation
                if (is_exp) printf("%6.2f %6.2f %6.2f %7.5f\n",
                                   xi+ri_slice*cos(c),yi+ri_slice*sin(c),z_slice,c);
            }
            printf("\n");
        }
#endif /* Debug */

    }
    return sasa;
}

static double sasa_sum_angles(int n_buried, double *restrict a, double *restrict b)
{
    /* Innermost function in L&R, could potentially be sped up, but
       probably requires rethinking, algorithmically. Perhaps
       recursion could be rolled out somehow. */
    char excluded[n_buried], n_exc = 0, n_overlap = 0;
    double ai,aj,bi,bj,d,inf,sup,buried_angle;
    memset(excluded,0,n_buried);

    for (int i = 0; i < n_buried; ++i) {
        if (excluded[i]) continue;
        for (int j = 0; j < n_buried; ++j) {
            if (i == j) continue;
            if (excluded[j]) continue;
            
            //check for overlap
            bi = b[i]; ai = a[i]; //will be updating throughout the loop
            bj = b[j], aj = a[j];
            for (;;) {
                d = bj - bi;
                if (d > M_PI) bj -= 2*M_PI;
                else if (d < -M_PI) bj += 2*M_PI;
                else break;
            }
            if (fabs(d) > ai+aj) continue;
            ++n_overlap;

            //calculate new joint interval
            inf = fmin(bi-ai,bj-aj);
            sup = fmax(bi+ai,bj+aj);
            b[i] = (inf + sup)/2.0;
            a[i] = (sup - inf)/2.0;
            if (a[i] > M_PI) return 0;
            if (b[i] > M_PI) b[i] -= 2*M_PI;
            if (b[i] < -M_PI) b[i] += 2*M_PI;

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
    buried_angle = 0;
    for (int i = 0; i < n_buried; ++i) {
        buried_angle += 2.0*a[i];
    }
    return 2*M_PI - buried_angle;
}

//discretizes circle, but for now it's not faster than other algorithm..
static double sasa_sum_angles_int(int n_buried, double *restrict a, double *restrict b)
{
    const int res = 360; // "degrees per circle", can be increased to improve accuracy
    char buried[res];
    const double rad2uni = 1./(2*M_PI);
    
    memset(buried,0,res);

    for (int i = 0; i < n_buried; ++i) {
        const double bi = b[i] + 2*M_PI; //makes sure there are no negative values
        const double ai = a[i];

        int inf = ((int) ((bi - ai) * rad2uni*res)) % res;
        int sup = ((int) ((bi + ai) * rad2uni*res)) % res;
        
        if (inf < sup) {
            memset(buried+inf,1,sup-inf);
        } else {
            memset(buried,1,sup);
            memset(buried+inf,1,res-inf);
        }
    }
    int count = 0;
    for (int i = 0; i < res; ++i) {
        if (buried[i]) ++count;
    }
    return 2*M_PI*((res - count)/(double)(res));
}

