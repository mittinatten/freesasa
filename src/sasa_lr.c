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

#ifdef __GNUC__
#define __attrib_pure__ __attribute__((pure))
#else
#define __attrib_pure__
#endif

extern const char *freesasa_name;
extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

//calculation parameters (results stored in *sasa)
typedef struct {
    int n_atoms;
    const double *radii;
    const freesasa_coord_t *xyz;
    const int **nb; // neighbors
    const int *nn; // number of neighbors
    const double **nb_xyd; // neighbours, xy-distance
    const double **nb_xd; // neighbours, x-distance
    const double **nb_yd; // neighbours, y-distance
    double delta; // slice width
    double min_z; // bounds of the molecule
    double max_z;
    double *sasa; // results
} sasa_lr_t; 

typedef struct {
    int n_slice; //number of atoms in slice
    double z; //the mid-point of the slice
    char * restrict in_slice;
    int * restrict idx; //index in slice to global numbering
    int * restrict xdi; //global numbering to index in slice
    double * restrict DR; //corrective multiplicative factor (D in L&R paper)
    double * restrict r; //radius in slice;
} sasa_lr_slice_t;

#if HAVE_LIBPTHREAD
static void sasa_lr_do_threads(int n_threads, sasa_lr_t*);
static void *sasa_lr_thread(void *arg);
#endif

static void sasa_add_slice_area(double z, sasa_lr_t*);

// the z argument is only really necessary for the debugging section
static void sasa_exposed_arcs(sasa_lr_slice_t *,
			      double * restrict exposed_arc,
			      const sasa_lr_t *);

/** a and b are a set of alpha and betas (in the notation of the
    manual). This function finds the union of those intervals on the
    circle, and returns 2*PI minus the length of the joined
    interval(s) (i.e. the exposed arc length). Does not necessarily
    leave a and b in a consistent state. */
static double sasa_sum_angles(int n_buried, double *a, double *b);

/** Calculate contacts, given coordinates and radii. The array nb will
    have a list of neighbors to each atom, nn will say how many
    neighbors each atom has, and nb_xyd the distance between
    neighbours in the xy-plane. The arrays nn, nb and nb_xyd should be
    of size n_atoms. The elements of nb and nb_xyd are dynamically
    allocated to be of size nn[i]. **/
static void sasa_get_contacts(int **nb, int *nn, double **nb_xyd,
			      double **nb_xd, double **nb_yd,
                              const freesasa_coord_t *xyz, const double *radii);


int freesasa_lee_richards(double *sasa,
			 const freesasa_coord_t *xyz,
			 const double *atom_radii,
			 double probe_radius,
			 double delta,
			 int n_threads)
{
    /* Steps:
       Define slice range
       For each slice:
       1. Identify member atoms
       2. Calculate their radii in slice
       3. Calculate exposed arc-lengths for each atom
       Sum up arc-length*delta for each atom
    */
    size_t n_atoms = freesasa_coord_n(xyz);
    int return_value = FREESASA_SUCCESS;
    if (n_atoms == 0) {
	return freesasa_warn("Attempting Lee & Richards calculation "
			    "on empty coordinates");
    }
    // determine slice range and init radii and sasa arrays
    double max_z=-1e50, min_z=1e50;
    double max_r = 0;
    double radii[n_atoms];
    const double *v = freesasa_coord_all(xyz);
    for (size_t i = 0; i < n_atoms; ++i) {
        radii[i] = atom_radii[i] + probe_radius;
        double z = v[3*i+2], r = radii[i];
        max_z = z > max_z ? z : max_z;
        min_z = z < min_z ? z : min_z;
        sasa[i] = 0;
        max_r = r > max_r ? r : max_r;
    }
    min_z -= max_r;
    max_z += max_r;
    min_z += 0.5*delta;
 
    // determine which atoms are neighbours
    int *nb[n_atoms], nn[n_atoms];
    double *nb_xyd[n_atoms], *nb_xd[n_atoms], *nb_yd[n_atoms];
    sasa_get_contacts((int**)nb, (int*)nn, (double**)nb_xyd, 
		      (double**)nb_xd, (double**)nb_yd, xyz, radii);
    sasa_lr_t lr = {.n_atoms = n_atoms, .radii = radii, .xyz = xyz,
                    .nb = (const int**)nb, .nn = nn, 
		    .nb_xyd = (const double**)nb_xyd, 
		    .nb_xd = (const double**)nb_xd, 
		    .nb_yd = (const double**)nb_yd, 
		    .delta = delta, 
                    .min_z = min_z, .max_z = max_z, .sasa = sasa};
    
    if (n_threads > 1) {
#if HAVE_LIBPTHREAD
        sasa_lr_do_threads(n_threads, &lr);
#else
        return_value = freesasa_warn("program compiled for single-threaded use, "
                                    "but multiple threads were requested. Will "
                                    "proceed in single-threaded mode.\n");
        n_threads = 1;
#endif
    } 
    if (n_threads == 1) {
        // loop over slices
        for (double z = min_z; z < max_z; z += delta) {
            sasa_add_slice_area(z,&lr);
        }
    }    
    for (int i = 0; i < n_atoms; ++i) {
	free(nb[i]);
	free(nb_xyd[i]);
	free(nb_xd[i]);
	free(nb_yd[i]);
    }
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
	} else {
            slice->in_slice[i] = 0;
	    slice->xdi[i] = -1;
        }
    }
    slice->n_slice = n_slice;
    return slice;
}

static void sasa_free_slice(sasa_lr_slice_t* slice)
{
    free(slice->idx);
    free(slice->xdi);
    free(slice->in_slice);
    free(slice->r);
    free(slice->DR);
    free(slice);
}

static void sasa_add_slice_area(double z, sasa_lr_t *lr)
{
    sasa_lr_slice_t *slice = sasa_init_slice(z,lr);
    double * restrict exposed_arc = (double* )malloc(slice->n_slice*(sizeof(double)));
    
    //find exposed arcs
    sasa_exposed_arcs(slice, exposed_arc, lr);
    
    // calculate contribution to each atom's SASA from the present slice
    for (int i = 0; i < slice->n_slice; ++i) {
        lr->sasa[slice->idx[i]] += exposed_arc[i]*slice->r[i]*slice->DR[i];
    }

    sasa_free_slice(slice);
}

static void sasa_exposed_arcs(sasa_lr_slice_t *slice, 
			      double * restrict exposed_arc, 
			      const sasa_lr_t *lr)
{
    const int n_slice = slice->n_slice;
    const int * restrict nn = lr->nn;
    const int ** restrict nb = lr->nb;
    const double ** restrict nb_xyd = lr->nb_xyd;
    const double * restrict r = slice->r;
    
    char is_completely_buried[n_slice]; // keep track of completely buried circles
    for (int i = 0; i < n_slice; ++i) is_completely_buried[i] = 0;

    //loop over atoms in slice
    //lower-case i,j is atoms in the slice, upper-case I,J are their
    //corresponding global indexes
    for (int i = 0; i < n_slice; ++i) {
        if (is_completely_buried[i]) {
            continue;
        }
        int I = slice->idx[i];
        double ri = slice->r[i], a[n_slice], b[n_slice];
        int n_buried = 0;
	exposed_arc[i] = 0;
        // loop over neighbors 
        for (int ni = 0; ni < nn[I]; ++ni) {
            int J = nb[I][ni];
	    if (slice->in_slice[J] == 0) continue;
            assert (I != J);
	    int j = slice->xdi[J];
            double rj = r[j]; //, xij = v[3*J]-xi, yij = v[3*J+1]-yi;
            double d = nb_xyd[I][ni]; //sqrt(xij*xij+yij*yij);
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
        if (is_completely_buried[i] == 0) 
            exposed_arc[i] = sasa_sum_angles(n_buried,a,b);
        
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

static double sasa_sum_angles(int n_buried, double *a, double *b)
{
    /* Innermost function in L&R, could potentially be sped up, but
       probably requires rethinking, algorithmically. Perhaps
       recursion could be rolled out somehow. */
    int excluded[n_buried], n_exc = 0, n_overlap = 0;
    for (int i = 0; i < n_buried; ++i)  {
        excluded[i] = 0;
        assert(a[i] > 0);
    }
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

static void sasa_get_contacts(int **nb, int *nn, double **nb_xyd, 
			      double **nb_xd, double **nb_yd,
                              const freesasa_coord_t *xyz, const double *radii)
{
    /* For low resolution L&R this function is the bottleneck in
       speed. Will also depend on number of atoms. */
    size_t n_atoms = freesasa_coord_n(xyz);
    for (int i = 0; i < n_atoms; ++i) {
        nn[i] = 0;
        nb[i] = NULL;
	nb_xyd[i] = NULL;
	nb_xd[i] = NULL;
	nb_yd[i] = NULL;
    }
    const double *restrict v = freesasa_coord_all(xyz);

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
}

