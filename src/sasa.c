#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "sasa.h"

#include "shrake_rupley_points.h"

void sasa_exposed_angles(int n_slice, double *x, double *y, double z, double *r, double *exposed_angle, int **nb, int *nn);
double sasa_sum_angles(int n_buried, double *a, double *b);

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
        //paranthesis too make sure floating point division is used (I am paranoid about this).
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
    double max_z=-1e50, min_z=1e50;
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

    //allocate some arrays to keep track of things
    int **contact = (int**) malloc(sizeof(int*)*n_atoms);
    int **nb = (int**) malloc(sizeof(int*)*n_atoms);
    for (int i = 0; i < n_atoms; ++i) {
        contact[i] = (int*) malloc(sizeof(int)*n_atoms);
        nb[i] = (int*) malloc(sizeof(int)*200); //assume atom won't have more than 200 neighbours;
    }

    // determine which atoms are neighbours (speeds up calculations a bit later on (factor 2 or so)).
    for (int i = 0; i < n_atoms; ++i) {
        double ri = radii[i]+PROBE_RADIUS;
        for (int j = i+1; j < n_atoms; ++j) {
            double rj = radii[j] + PROBE_RADIUS;
            double cut = ri + rj;
            if (vector3_dist2(&xyz[i],&xyz[j]) < cut*cut) {
                contact[i][j] = contact[j][i] = 1;
            } else {
                contact[i][j] = contact[j][i] = 0;
            }
        }
    }
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
        int nn[n_atoms];
        for (int i = 0; i < n_slice; ++i) nn[i] = 0;
        for (int i = 0; i < n_slice; ++i) {
            for (int j = i+1; j < n_slice; ++j) {
                if (contact[idx[i]][idx[j]]) {
                    nb[i][nn[i]++] = j;
                    nb[j][nn[j]++] = i;
                }
            }
            assert(nn[i] < 200 && "An atom had 200 or more neighbors, this should not be possible.");
        }
        sasa_exposed_angles(n_slice, x, y, z, r, exposed_angle, nb, nn);
        // calculate contribution to each atom's SASA from the present slice
        for (int i = 0; i < n_slice; ++i) {
            sasa[idx[i]] += exposed_angle[i]*r[i]*delta;
            //printf("s[%d]: %f\n",idx[i],sasa[idx[i]]);
        }
    }
    for (int i = 0; i < n_atoms; ++i) {
        free(contact[i]);
        free(nb[i]);
    }
    free(contact);
    free(nb);
}

void sasa_exposed_angles(int n_slice, double *x, double *y, double z, double *r,
                         double *exposed_angle, int **nb, int *nn)
{
    int is_buried[n_slice]; // keep track of completely buried circles
    for (int i = 0; i < n_slice; ++i) is_buried[i] = 0;
    //loop over atoms in slice
    for (int i = 0; i < n_slice; ++i) {
        double ri = r[i], a[n_slice], b[n_slice];
        int n_buried = 0;
        exposed_angle[i] = 0;
        if (is_buried[i]) {
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
            if (d + ri < rj) {is_buried[i] = 1; break;} // circle i is completely inside j
            if (d + rj < ri) {is_buried[j] = 1; continue;} // circle j is completely inside i

            // half the arclength occluded from circle i due to verlap with circle j
            double alpha = acos ((ri*ri + d*d - rj*rj)/(2.0*ri*d));
            // the polar coordinates angle of the vector connecting i and j
            double beta = atan2 (yij,xij);

            a[n_buried] = alpha;
            b[n_buried] = beta;

            ++n_buried;
        }
	if (is_buried[i] != 0) exposed_angle[i] = sasa_sum_angles(n_buried,a,b);

#ifdef DEBUG
        if (is_buried != 0) {
            for (double c = 0; c < 2*PI; c += PI/45.0) {
                int is_exp = 1;
                for (int i = 0; i < n_buried; ++i) {
                    if ((c > b[i]-a[i] && c < b[i]+a[i]) ||
                        (c - 2*PI > b[i]-a[i] && c - 2*PI < b[i]+a[i]) ||
                        (c + 2*PI > b[i]-a[i] && c + 2*PI < b[i]+a[i])) { is_exp = 0; break; }
                }
                // print the arcs used in calculation
                if (is_exp) printf("%6.2f %6.2f %6.2f %7.5f\n",x[i]+ri*cos(c),y[i]+ri*sin(c),z,c);
		//if (is_exp) exposed_angle[i] += PI/45.0;
	    }
        }
#endif
    }
}

double sasa_sum_angles(int n_buried, double *a, double *b)
{
    int excluded[n_buried];
    //static int i = 0; ++i; printf("\n%d\n",i);
    for (int i = 0; i < n_buried; ++i)  {
	excluded[i] = 0;
	assert(a[i] > 0);
	//printf("[%4.2f, %4.2f] ",(b[i]-a[i])/PI,(b[i]+a[i])/PI);
    }
    //printf("\n");
   
    for (int i = 0; i < n_buried; ++i) {
        if (excluded[i]) continue;
        for (int j = 0; j < n_buried; ++j) {
	    if (excluded[j]) continue;
	    if (i == j) continue;
            double d;
	    double bi = b[i], ai = a[i]; //will be updating throughout the loop
	    double bj = b[j], aj = a[j];
	    //printf ("{%4.2f, %4.2f, %4.2f, %4.2f}\n",(bi-ai)/PI,(bi+ai)/PI,(bj-aj)/PI,(bj+aj)/PI);
            for(;;) {
                d = bj - bi;
                if (d > PI) bj -= 2*PI;
                else if (d < -PI) bj += 2*PI;
                else break;
            }
            if (fabs(d) > ai+aj) continue;
            double inf_i = bi-ai, inf_j = bj-aj;
            double sup_i = bi+ai, sup_j = bj+aj;
            double inf = inf_i < inf_j ? inf_i : inf_j;
            double sup = sup_i > sup_j ? sup_i : sup_j;
	    //printf("%4.2f %4.2f %4.2f\n", d/PI, inf/PI, sup/PI);
            b[i] = (inf + sup)/2.0;
            a[i] = fabs(sup - inf)/2.0;
            b[j] = a[j] = 0;
	    if (a[i] > PI) return 0;
            if (b[i] > PI) b[i] -= 2*PI;
            if (b[i] < -PI) b[i] += 2*PI;
            excluded[j] = 1;
        }
    }
    double buried_angle = 0;
    for (int i = 0; i < n_buried; ++i) {
	//if (excluded[i] == 0) // should not be used!! 
	buried_angle += 2.0*a[i];
	//if (a[i] > 0 ) 
	//    printf("[%4.2f, %4.2f] ",(b[i]-a[i])/PI,(b[i]+a[i])/PI);
    }
    double exposed_angle = 2*PI - buried_angle;
    if (exposed_angle < 0) exposed_angle = 0; // will be the case for almost all buried atoms
    
    //printf("%4.2f\n",exposed_angle/PI);
    //if (i > 992) exit(0);

    return exposed_angle;
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
