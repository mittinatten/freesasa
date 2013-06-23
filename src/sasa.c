#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "sasa.h"

#include "shrake_rupley_points.h"

void sasa_exposed_angles(int n_slice, double *x, double *y, double *r, double *exposed_angle);
double sasa_calc_exposed_angle(int n_buried, double *a, double *b);

void sasa_shrake_rupley(double *sasa,
                        const vector3 *xyz,
                        const double *r,
                        size_t n_atoms,
                        int n_points)
{
    assert(n_points <= MAX_SR_POINTS);
    const double *srp = shrake_rupley_points;

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
        //paranthesis too make sure floating point division is used.
        sasa[i] = (4.0*PI*ri*ri*n_surface)/n_points;
    }
}

void sasa_lee_richards(double *sasa,
                       const vector3 *xyz,
                       const double *radii,
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
    double max_z=-1e50, min_z=1e50, start_z;
    double max_r = 0;
    for (size_t i = 0; i < n_atoms; ++i) {
	double z = xyz[i].z, r = radii[i];
	max_z = z > max_z ? z : max_z;
	min_z = z < min_z ? z : min_z;
	sasa[i] = 0;
	max_r = r > max_r ? r : max_r;
    }
    min_z -= max_r;
    max_z += max_r;

    // loop over slices
    for (double z = min_z + 0.5*delta; z < max_z; z += delta) {
	double x[n_atoms], y[n_atoms], r[n_atoms];
	int n_slice = 0;
	double exposed_angle[n_atoms];
	int idx[n_atoms];
	// locate atoms in each slice
	for (size_t i = 0; i < n_atoms; ++i) {
	    double ri = radii[i] + PROBE_RADIUS;
	    double d = fabs(xyz[i].z-z);
	    if (d < ri) {
		x[n_slice] = xyz[i].x;
		y[n_slice] = xyz[i].y;
		r[n_slice] = sqrt(ri*ri-d*d);
		idx[n_slice] = i;
		++n_slice;
	    }
	}
	sasa_exposed_angles(n_slice, x, y, r, exposed_angle);
	// calculate contribution to each atom's SASA from the present slice
	for (int i = 0; i < n_slice; ++i) {
	    sasa[idx[i]] += exposed_angle[i]*r[i]*delta; 
	    //printf("s[%d]: %f\n",idx[i],sasa[idx[i]]);
	}
    }    
}

void sasa_exposed_angles(int n_slice, double *x, double *y, double *r, double *exposed_angle)
{
    for (int i = 0; i < n_slice; ++i) {
	double ri = r[i];
	double a[n_slice], b[n_slice];
	int n_buried = 0;
	// loop over atoms in slice
	for (int j = 0; j < n_slice; ++j) {
	    if (i == j) continue;
	    double rj = r[j];
	    double xij = x[j]-x[i];
	    double yij = y[j]-y[i];
	    double rij = sqrt(xij*xij+yij*yij);
	    
	    // reasons to skip calculation
	    if (rij > ri + rj) continue;     // atoms aren't in contact
	    if (rij+ri < rj) break; // circle i is completely inside j
	    if (rij+rj < ri) continue; // circle j is completely inside i
	    
	    // half the arclength occluded from circle i
	    double alpha = acos ((ri*ri+rij*rij-rj*rj)/(2.0*ri*rij));
	    // the polar coordinates angle of the vector connecting i and j
	    double gamma = atan2 (yij,xij);

	    a[n_buried] = gamma - alpha;
	    b[n_buried] = gamma + alpha; 
	    ++n_buried;
	}
	exposed_angle[i] = sasa_calc_exposed_angle(n_buried,a,b);
    }
}

// does not leave a and b arrays unchanged
double sasa_calc_exposed_angle(int n_buried, double *a, double *b)
{
    // loop over buried arcs
    for (int j = 0; j < n_buried; ++j) {
	double bj = b[j] > 2*PI ? b[j]-2*PI : b[j];
	for (int k = 0; k < n_buried; ++k) {
	    if (j == k) continue;
	    double bk = b[k] > 2*PI ? b[k]-2*PI : b[k];
	    // overlap condition (multiples of 2*PI already handled)
	    if (a[j] < bk && bj > a[k]) {
		// interval j is the union between j and k, and k is made empty.
		a[j] = a[j] < a[k] ? a[j] : a[k]; 
		b[j] = b[j] > b[k] ? b[j] : b[k]; // important not to use bj and bk here
		a[k] = b[k] = 0;
	    }
	}
    }
    double arc = 0;
    for (int j = 0; j < n_buried; ++j) arc += b[j]-a[j];
    arc = 2*PI - arc; // the exposed arc
    if (arc < 0) arc = 0; // will be the case for almost all buried atoms
    return arc;
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
    for (int i = 0; i < nc; ++i) {
        fprintf(out,"%s\t%6.2f\n",ac.class2str[i],s[i]);
    }
}
