/*
  Copyright Simon Mitternacht 2013.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sasa.h"
#include "srp.h"

/** internal functions for L&R calculations **/
void sasa_add_slice_area(double *sasa, double z, double delta, const vector3 *xyz, 
			 const double *radii, const int **nb, const int *nn, 
			 int n_atoms);
// the z argument is only really necessary for the debugging section
void sasa_exposed_arcs(int n_slice, const double *x, const double *y, double z, 
			 const double *r, double *exposed_arc, 
			 const int **nb, const int *nn);
//does not necessarily leave a and b in a consistent state
double sasa_sum_angles(int n_buried, double *a, double *b);

/** Calculate contacts, given coordinates and radii. The array nb will
    have a list of neighbors to each atom, nn will say how many
    neighbors each atom has. The arrays nn and nb should be of size
    n_atoms. The elements of n_atoms are dynamically allocated to be
    of size nn[i]. **/
void sasa_get_contacts(int **nb, int *nn, int n_atoms, 
		       const vector3 *xyz, const double *radii);

void sasa_shrake_rupley(double *sasa,
                        const vector3 *xyz,
                        const double *r,
                        size_t n_atoms,
                        int n_points)
{
    const double *srp = srp_get_points(n_points);

    //turn test-points into vector3's
    vector3 srp_v3[n_points];
    for (int k = 0; k < n_points; ++k) {
        vector3_set(&srp_v3[k],srp[k*3],srp[k*3+1],srp[k*3+2]);
        vector3_setlength(&srp_v3[k],1.0); //just to be sure
    }
    vector3 test;

    //calculate SASA
    for (int i = 0; i < n_atoms; ++i) {

        /* this array keeps track of which testpoints belonging to
           a certain atom overlap with other atoms */
        char spcount[n_points];
        for (int k = 0; k < n_points; ++k) {
            spcount[k] = 0;
        }

        double ri = r[i]+PROBE_RADIUS;
        for (int j = 0; j < n_atoms; ++j) {
            if (j == i) continue;
            const double rj = r[j]+PROBE_RADIUS;
            const double cut2 = (ri+rj)*(ri+rj);
            if (vector3_dist2(&xyz[i],&xyz[j]) > cut2)
                continue;
            for (int k = 0; k < n_points; ++k) {
                if (spcount[k]) continue;

                /* this probably cannot be done beforehand, because
                   changing the length repeatedly could cause
                   round-off errors */
                vector3_copy(&test, &srp_v3[k]);

                // test-point vectors have length 1.0
                vector3_multiply(&test, ri);
                vector3_add(&test, &xyz[i]);

                if (vector3_dist2(&test, &xyz[j]) <= rj*rj) {
                    spcount[k] = 1;
                }
                /* i.e. if |xyz[i]+ri*srp[k] - xyz[j]| <= rj we have an
                   overlap. */
            }
        }
        int n_surface = 0;
        for (int k = 0; k < n_points; ++k) {
            if (!spcount[k]) ++n_surface;
#ifdef DEBUG
            if (!spcount[k]) {
                vector3_copy(&test,&srp_v3[k]);
                vector3_multiply(&test,ri);
                vector3_add(&test,&xyz[i]);
                printf("%f %f %f\n",test.x,test.y,test.z);
            }
#endif
        }
        //paranthesis to make sure floating point division is used (I am paranoid about this).
        sasa[i] = (4.0*PI*ri*ri*n_surface)/n_points;
    }
}

void sasa_lee_richards(double *sasa,
                       const vector3 *xyz,
                       const double *atom_radii,
                       size_t n_atoms,
                       double delta)
{
    /* Steps:
       Define slice range
       For each slice:
       1. Identify member atoms
       2. Calculate their radii in slice
       3. Calculate exposed arc-lengths for each atom
       Sum up arc-length*delta for each atom
    */

    // determine slice range and init radii and sasa arrays
    double max_z=-1e50, min_z=1e50;
    double max_r = 0;
    double radii[n_atoms];
    for (size_t i = 0; i < n_atoms; ++i) {
	radii[i] = atom_radii[i] + PROBE_RADIUS;
        double z = xyz[i].z, r = radii[i];
        max_z = z > max_z ? z : max_z;
        min_z = z < min_z ? z : min_z;
        sasa[i] = 0;
        max_r = r > max_r ? r : max_r;
    }
    min_z -= max_r;
    max_z += max_r;

    
    // determine which atoms are neighbours
    int *nb[n_atoms], nn[n_atoms];
    sasa_get_contacts((int**)nb,(int*)nn,n_atoms,xyz,radii);
    
    // loop over slices
    for (double z = min_z + 0.5*delta; z < max_z; z += delta) {
	sasa_add_slice_area(sasa,z,delta,xyz,radii,(const int**)nb,(int*)nn,n_atoms);
	//could be trivially parallelized if summation is moved outside of above function
    }
    for (int i = 0; i < n_atoms; ++i) {
	free(nb[i]);
    }
}

void sasa_add_slice_area(double *sasa, double z, double delta, const vector3 *xyz, 
			 const double *radii, const int **nb, const int *nn, int n_atoms)
{
    double x[n_atoms], y[n_atoms], r[n_atoms], DR[n_atoms];
    int n_slice = 0;
    double exposed_arc[n_atoms];
    int idx[n_atoms], xdi[n_atoms], in_slice[n_atoms], nn_slice[n_atoms], *nb_slice[n_atoms];

    // locate atoms in each slice and do some initialization
    for (size_t i = 0; i < n_atoms; ++i) {
	double ri = radii[i];
	double d = fabs(xyz[i].z-z);
	if (d < ri) {
	    x[n_slice] = xyz[i].x; y[n_slice] = xyz[i].y;
	    r[n_slice] = sqrt(ri*ri-d*d);
	    //multiplicative factor when arcs are summed up later (according to L&R paper)
	    DR[n_slice] = ri/r[n_slice]*(delta/2. +
			      (delta/2. < ri-d ? delta/2. : ri-d));
	    idx[n_slice] = i;
	    xdi[i] = n_slice;
	    ++n_slice;
	    in_slice[i] = 1;
	} else {
	    in_slice[i] = 0;
	}
    }
    for (int i = 0; i < n_slice; ++i) { 
	nn_slice[i] = 0;
	exposed_arc[i] = 0;
    }

    for (int i = 0; i < n_slice; ++i) {
	int i2 = idx[i];
	int j2;
	for (int j = 0; j < nn[i2]; ++j) {
	    if (in_slice[j2 = nb[i2][j]]) {
		++nn_slice[i];
		nb_slice[i] = realloc(nb_slice[i],sizeof(int)*nn_slice[i]);
		nb_slice[i][nn_slice[i]-1] = xdi[j2];
	    }
	}
    }

    //find exposed arcs
    sasa_exposed_arcs(n_slice, x, y, z, r, exposed_arc, (const int**)nb_slice, nn_slice);
    
    // calculate contribution to each atom's SASA from the present slice
    for (int i = 0; i < n_slice; ++i) {
	sasa[idx[i]] += exposed_arc[i]*r[i]*DR[i];
    }
}
void sasa_exposed_arcs(int n_slice, const double *x, const double *y, double z, const double *r,
		       double *exposed_arc, const int **nb, const int *nn)
{
    int is_completely_buried[n_slice]; // keep track of completely buried circles
    for (int i = 0; i < n_slice; ++i) is_completely_buried[i] = 0;
    //loop over atoms in slice
    for (int i = 0; i < n_slice; ++i) {
        double ri = r[i], a[n_slice], b[n_slice];
        int n_buried = 0;
        exposed_arc[i] = 0;
        if (is_completely_buried[i]) {
            continue;
        }
        // loop over neighbors in slice
        for (int ni = 0; ni < nn[i]; ++ni) {
	    int j = nb[i][ni];
            assert (i != j);
	    double rj = r[j], xij = x[j]-x[i], yij = y[j]-y[i];
            double d = sqrt(xij*xij+yij*yij);
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
            double beta = atan2 (yij,xij);

            a[n_buried] = alpha;
            b[n_buried] = beta;

            ++n_buried;
        }
	if (is_completely_buried[i] == 0) 
	    exposed_arc[i] = sasa_sum_angles(n_buried,a,b);

#ifdef DEBUG
	if (is_completely_buried[i] == 0) {
	    //exposed_arc[i] = 0;
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
				   x[i]+ri*cos(c),y[i]+ri*sin(c),z,c);
		//if (is_exp) exposed_arc[i] += PI/45.0;
	    }
	    printf("\n");
        }
#endif
    }
}

double sasa_sum_angles(int n_buried, double *a, double *b)
{
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

void sasa_per_atomclass(FILE *out, atomclassifier ac,
                        protein *p, double *sasa)
{
    int nc = ac.nclasses;
    double s[nc];
    for (int i = 0; i < nc; ++i) {
        s[i] = 0;
    }
    for (size_t i = 0; i < protein_n(p); ++i) {
        int class = ac.classify(protein_atom_res_name(p,i),
                                protein_atom_name(p,i));
        s[class] += sasa[i];
    }
    fprintf(out,"\n> %s\n",ac.name);
    for (int i = 0; i < nc; ++i) {
	if (s[i] > 0.0) fprintf(out,"%s %9.2f\n",ac.class2str[i],s[i]);
    }
}

void sasa_get_contacts(int **nb, int *nn, int n_atoms, 
		       const vector3 *xyz, const double *radii)
{
    /* For low resolution L&R this function is the bottleneck in
       speed. Will also depend on number of atoms. */
    
    for (int i = 0; i < n_atoms; ++i) {
	nn[i] = 0;
	nb[i] = NULL;
    }

    for (int i = 0; i < n_atoms; ++i) {
        double ri = radii[i];
        const vector3 *vi = &xyz[i];
	//double xi = vi->x, yi = vi->y, zi = vi->z;
	for (int j = i+1; j < n_atoms; ++j) {
            double rj = radii[j];
	    double cut2 = (ri+rj)*(ri+rj);
	    const vector3 *vj = &xyz[j];
	    
	    /* could be minor speed improvement, most pairs of atoms will be
	       far away from each other on at least one axis */
	    /*
	    double xj = vj->x, yj = vj->y, zj = vj->z;
            if ((xj-xi)*(xj-xi) > cut2 ||
		(yj-yi)*(yj-yi) > cut2 ||
		(zj-zi)*(zj-zi) > cut2) {
		continue;
		}*/
	    if (vector3_dist2(vi,vj) < cut2) {
		++nn[i]; ++nn[j];
		nb[i] = realloc(nb[i],sizeof(int)*nn[i]);
		nb[j] = realloc(nb[j],sizeof(int)*nn[j]);
		nb[i][nn[i]-1] = j;
		nb[j][nn[j]-1] = i;
	    }
	}
    }

}

