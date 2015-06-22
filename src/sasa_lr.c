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
    int n_slice; //number of atoms in slice
    double z; //the mid-point of the slice
    char *in_slice;
    int *idx; //index in slice to global numbering
    int *xdi; //global numbering to index in slice
    double *DR; //corrective multiplicative factor (D in L&R paper)
    double *r; //radius in slice;
    double *exposed_arc; //exposed arc length (in Å) for each atom
} sasa_lr_slice;

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr*);
static void *sasa_lr_thread(void *arg);
#endif

/** Init slice and sum up arc-length */
static void sasa_add_slice_area(double z, sasa_lr*);

/** Find which arcs are exposed in a slice */
static void sasa_exposed_arcs(sasa_lr_slice * restrict,
                              const sasa_lr *);

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

static void sasa_add_atoms(sasa_lr *lr);

/** Initialize object to be used for L&R calculation */
static sasa_lr* freesasa_init_lr(double *sasa,
                                   const freesasa_coord *xyz,
                                   const double *atom_radii,
                                   double probe_radius,
                                   double delta) {
    sasa_lr* lr = malloc(sizeof(sasa_lr));
    assert(lr);
    const int n_atoms = freesasa_coord_n(xyz);
    double max_z=-1e50, min_z=1e50;
    double max_r = 0;
    double *radii = malloc(sizeof(double)*n_atoms);
    assert(radii);
    const double *v = freesasa_coord_all(xyz);
    //find bounds of protein along z-axis and init radii
    for (int i = 0; i < n_atoms; ++i) {
        radii[i] = atom_radii[i] + probe_radius;
        double z = v[3*i+2], r = radii[i];
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

    if (freesasa_coord_n(xyz) == 0) {
        return freesasa_warn("%s: empty coordinates",__func__);
    }

    // determine slice range and init radii and sasa arrays,
    // and find which atoms are neighbors
    sasa_lr *lr = freesasa_init_lr(sasa, xyz, atom_radii,
                                   probe_radius, delta);

    // determine which atoms are neighbours
    lr->adj = freesasa_adjacency_new(xyz,lr->radii);

    //sasa_get_contacts(lr);

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
        // loop over slices
        sasa_add_atoms(lr);
        
        /*for (double z = lr->min_z; z < lr->max_z; z += lr->delta) {
            sasa_add_slice_area(z,lr);
            }*/
    }
    freesasa_free_lr(lr);

    return return_value;
}

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr *lr)
{
    double *t_sasa[n_threads];
    pthread_t thread[n_threads];
    sasa_lr lrt[n_threads];
    const double max_z = lr->max_z, min_z = lr->min_z, delta = lr->delta;
    int n_slices = (int)ceil((max_z-min_z)/delta);
    int n_perthread = n_slices/n_threads;
    for (int t = 0; t < n_threads; ++t) {
        t_sasa[t] = malloc(sizeof(double)*lr->n_atoms);
        assert(t_sasa[t]);
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
    sasa_lr *lr = ((sasa_lr*) arg);
    for (double z = lr->min_z; z < lr->max_z; z += lr->delta) {
        sasa_add_slice_area(z, lr);
    }
    pthread_exit(NULL);
}
#endif /* pthread */

// find which atoms are in slice
static sasa_lr_slice* sasa_init_slice(double z, const sasa_lr *lr)
{
    sasa_lr_slice *slice = malloc(sizeof(sasa_lr_slice));
    assert(slice);
    int n_atoms = lr->n_atoms;
    int chunk = n_atoms/10 > 16 ? n_atoms/10 : 16;
    int capacity = chunk;
    int n_slice = slice->n_slice = 0;
    double delta = lr->delta;
    double d_half = delta/2.;
    const double * restrict v = freesasa_coord_all(lr->xyz);
    const double * restrict r = lr->radii;

    slice->xdi = malloc(n_atoms*sizeof(int));
    assert(slice->xdi);
    slice->in_slice = malloc(n_atoms*sizeof(int));
    assert(slice->in_slice);
    slice->idx = malloc(chunk*sizeof(int));
    slice->r = malloc(chunk*sizeof(double));
    slice->DR = malloc(chunk*sizeof(double));
    assert(slice->idx); assert(slice->r); assert(slice->DR);
    slice->z = z;

    memset(slice->in_slice,0,sizeof(int)*n_atoms);

    // locate atoms in each slice and do some initialization
    for (int i = 0; i < n_atoms; ++i) {
        double ri = r[i];
        double d = fabs(v[3*i+2]-z);
        double r;
        if (d < ri) {
            if (n_slice+1 > capacity) {
                capacity += chunk;
                slice->idx = (int*) realloc(slice->idx,(capacity)*sizeof(int));
                slice->r = (double*) realloc(slice->r,(capacity)*sizeof(double));
                slice->DR = (double*) realloc(slice->DR,(capacity)*sizeof(double));
                assert(slice->idx && slice->r && slice->DR);
            }
            
            slice->idx[n_slice] = i;
            slice->xdi[i] = n_slice;
            slice->in_slice[i] = 1;
            slice->r[n_slice] = r = sqrt(ri*ri-d*d);
            slice->DR[n_slice] = ri/r*(d_half +
                                       (d_half < ri-d ? d_half : ri-d));
            ++n_slice;
        }  else {
            slice->xdi[i] = -1;
        } 
    }
    slice->n_slice = n_slice;
    slice->exposed_arc = malloc(n_slice*(sizeof(double)));
    assert(slice->exposed_arc);
    return slice;
}

static void sasa_free_slice(sasa_lr_slice* slice)
{
    if (slice) {
        free(slice->idx);
        free(slice->xdi);
        free(slice->in_slice);
        free(slice->r);
        free(slice->DR);
        free(slice->exposed_arc);
        free(slice);
    }
}

static void sasa_add_slice_area(double z, sasa_lr *lr)
{
    sasa_lr_slice *slice = sasa_init_slice(z,lr);

    //find exposed arcs
    sasa_exposed_arcs(slice, lr);

    // calculate contribution to each atom's SASA from the present slice
    for (int i = 0; i < slice->n_slice; ++i) {
        lr->sasa[slice->idx[i]] += slice->exposed_arc[i];
    }
    sasa_free_slice(slice);
}

static void sasa_exposed_arcs(sasa_lr_slice * restrict slice,
                              const sasa_lr *lr)
{
    const int n_slice = slice->n_slice;
    const freesasa_adjacency * restrict adj = lr->adj;
    const int *restrict nn = adj->nn;
    const double *restrict r = slice->r;
    const char *restrict in_slice = slice->in_slice;
    const int *restrict xdi = slice->xdi;
    char is_completely_buried[n_slice]; // keep track of completely buried circles
    memset(is_completely_buried,0,n_slice);
    
    //loop over atoms in slice
    //lower-case i,j is atoms in the slice, upper-case I,J are their
    //corresponding global indexes
    for (int i = 0; i < n_slice; ++i) {
        slice->exposed_arc[i] = 0;
        if (is_completely_buried[i]) {
            continue;
        }
        int I = slice->idx[i];
        double ri = r[i], a[nn[I]], b[nn[I]];
        int n_buried = 0;
        const int * restrict nbI = adj->nb[I];
        const double * restrict xydI = adj->nb_xyd[I];
        // loop over neighbors
        for (int ni = 0; ni < nn[I]; ++ni) {
            int J = nbI[ni];
            if (in_slice[J] == 0) continue;
            int j = xdi[J];
            double rj = r[j];
            double d = xydI[ni];
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
            double beta = atan2 (adj->nb_yd[I][ni],adj->nb_xd[I][ni]);
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
            for (double c = 0; c < 2*M_PI; c += M_PI/45.0) {
                int is_exp = 1;
                for (int i = 0; i < n_buried; ++i) {
                    if ((c > b[i]-a[i] && c < b[i]+a[i]) ||
                        (c - 2*M_PI > b[i]-a[i] && c - 2*M_PI < b[i]+a[i]) ||
                        (c + 2*M_PI > b[i]-a[i] && c + 2*M_PI < b[i]+a[i])) {
                        is_exp = 0; break;
                    }
                }
                // print the arcs used in calculation
                if (is_exp) printf("%6.2f %6.2f %6.2f %7.5f\n",
                                   xi+ri*cos(c),yi+ri*sin(c),slice->z,c);
            }
            printf("\n");
        }
#endif /* Debug */
    }
}

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
    int bottom = ((zi-ri)-lr->min_z)/delta;

    z0 = lr->min_z+bottom*delta;
    
    for (int j = 0; j < nni; ++j) {
        z_nb[j] = v[3*nbi[j]+2];
        r_nb[j] = r[nbi[j]];
    }
    
    for (z_slice = z0; z_slice < zi+ri; z_slice += delta) {
        double di = fabs(zi - z_slice);
        if (di > ri) continue; //deal with round off errors in z0
        double ri_slice = sqrt(ri*ri-di*di);
        double DR = ri/ri_slice*(d_half + fmin(d_half,ri-di));
        int n_buried = 0, completely_buried = 0;
        for (int j = 0; j < nni; ++j) {
            double zj = z_nb[j];
            double dj = fabs(zj - z_slice);
            double rj = r_nb[j];
            if (dj < rj) {
                double rj_slice = sqrt(rj*rj-dj*dj);
                //printf("di %f, ri %f, ri_s %f\n",dj,rj,rj_slice);
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

static void sasa_add_atoms(sasa_lr *lr)
{
    for (int i = 0; i < lr->n_atoms; ++i) {
        lr->sasa[i] = sasa_atom_area(lr,i); 
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
            if (i == j) continue;
            if (excluded[j]) continue;
            
            //check for overlap
            double bi = b[i], ai = a[i]; //will be updating throughout the loop
            double bj = b[j], aj = a[j];
            double d;
            for (;;) {
                d = bj - bi;
                if (d > M_PI) bj -= 2*M_PI;
                else if (d < -M_PI) bj += 2*M_PI;
                else break;
            }
            if (fabs(d) > ai+aj) continue;
            ++n_overlap;

            //calculate new joint interval
            double inf = fmin(bi-ai,bj-aj);
            double sup = fmax(bi+ai,bj+aj);
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
    double buried_angle = 0;
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
    memset(buried,0,res);

    const double rad2uni = 1./(2*M_PI);

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

