#include "sasa.h"

#include "shrake_rupley_points.c"

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
	vector3_set(&srp_v3[k],srp[i*3],srp[i*3+1],srp[i*3+2]);
    }
    
    //calculate SASA
    for (int i = 0; i < n_atoms; ++i) {
	
        /* this array keeps track of which testpoints belonging to
	   a certain atom overlap with other atoms */
	char spcount[n_points];
	for (int k = 0; k < n_points; ++k) {
	    spcount[k] = 0;
	}
	
	double ri = r[i];
	for (int j = 0; j < n_atoms; ++j) {
	    if (j == i) continue;
	    const double rj = r[j];
	    const double cut2 = (ri+rj)*(ri+rj);
	    if (vector3_dist2(xyz[i],xyz[j]) > cut2)
		continue;
	    vector3 test;
	    for (int k = 0; k < n_points; ++k) {
		if (spcount[k]) continue; 

                /* this probably cannot be done beforehand, because
                   changing the length repeatedly could cause
                   round-off errors */
		vector3_copy(&test, &spr_v3[k]); 

		vector3_setlength(&test, ri);
		vector3_add(&test, &xyz[i]);
		if (vector3_dist2(&test, &xyz[j]) <= rj*rj) {
		    spcount[k] = 1;
		}
		/* i.e. if |xyz[i]+ri*spr[k] - xyz[j]| <= rj we have an
		   overlap. */
	    }
	}
	int n_overlap = 0;
	for (int k = 0; k < n_points; ++k) {
	    if (spcount[k]) ++n_overlap;
	}
	sasa[i] = 4.0*PI*ri*ri*double(n_overlap)/double(n_points);
    }
}

double* sasa_lee_richards(double *sasa,
			  const vector3 *xyz,
			  const double *radii,
			  size_t n_atoms,
			  double grid)
{
    return NULL;
}
