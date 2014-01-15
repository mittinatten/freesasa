/*
  Copyright Simon Mitternacht 2013-2014.

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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "src/sasalib.h"
#include "src/coord.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

int number_tests;

double rel_err(double v1, double v2) {
    return fabs(v1-v2)/(fabs(v1)+fabs(v2));
}

void alg_failed(const char *alg,const char *test)
{
    printf("%s failed test: %s\n",alg,test);
}

double surface_hidden_sphere_intersection(double r1, double r2, double d) 
{
    if (d > r1 + r2) return 0;
    if (r1+d < r2) return 4*PI*r1*r1;
    if (r2+d < r1) return 4*PI*r2*r2;
    return PI/d * ( r1*(r2*r2 - (d-r1)*(d-r1))
		   +r2*(r1*r1 - (d-r2)*(d-r2)) ) ;
}
double surface_spheres_intersecting(double r1, double r2, double d)
{
    return 4*PI*(r1*r1+r2*r2) - surface_hidden_sphere_intersection(r1,r2,d);
}

double surface_two_spheres(const double *x, const double *r, double probe)
{
    double d2 = (x[0]-x[3])*(x[0]-x[3]) + (x[1]-x[4])*(x[1]-x[4]) 
	+ (x[2]-x[5])*(x[2]-x[5]);
    return surface_spheres_intersecting(r[0]+probe,r[1]+probe,sqrt(d2));
}

int test_sasa(sasalib_t *s, double ref, double tolerance,const char *test)
{
    sasalib_refresh(s);
    double err;
    if ((err = rel_err(ref,sasalib_area_total(s)))
	> tolerance) {
	alg_failed(sasalib_algorithm_name(s),test);
	printf("%f %f\n",ref,sasalib_area_total(s));
	return 1;
    }
    return 0;
}

int test_sasa_alg(sasalib_t *s, double tolerance) 
{
    int err = 0;
    int nt = 0;
    printf("Basic tests for algorithm '%s':\n",sasalib_algorithm_name(s));
    
    // Two spheres, compare with analytic results
    double coord[6] = {0,0,0,2,0,0};
    double r[2] = {1,2};
    double probe = sasalib_get_probe_radius(s);
    sasalib_link_coord(s,coord,r,2);
    
    ++nt;
    err += test_sasa(s,surface_two_spheres(coord,r,probe),tolerance,
		     "Two intersecting spheres along x-axis.");
    coord[3] = 0; coord[4] = 2;
    ++nt;
    err += test_sasa(s,surface_two_spheres(coord,r,probe),tolerance,
		     "Two intersecting spheres along y-axis."); 
    coord[4] = 0; coord[5] = 2;
    ++nt;
    err += test_sasa(s,surface_two_spheres(coord,r,probe),tolerance,
		     "Two intersecting spheres along z-axis.");

    // Four spheres in a plane, all calculations should give similar results
    double coord2[12] = {0,0,0, 1,0,0, 0,1,0, 1,1,0};
    double r2[4] = {1,1,2,1};
    sasalib_link_coord(s,coord2,r2,4);
    sasalib_refresh(s);
    double ref = sasalib_area_total(s);

    //translate
    for (int i = 0; i < 12; ++i) coord2[i] += 1.;
    ++nt;
    err += test_sasa(s,ref,tolerance,"Four spheres in plane, translated");

    //rotate 90 degrees round z-axis
    double coord3[12] = {0,1,0, 0,0,0, 1,1,0, 1,0,0};
    memcpy(coord2,coord3,12*sizeof(double));
    ++nt;
    err += test_sasa(s,ref,tolerance,"Four spheres in plane, rotated 90 deg round z-axis.");

    //rotate -45 degrees round z-axis
    double sqr2 = sqrt(2);
    double coord4[12] = {-1./sqr2,1./sqr2,0, 0,0,0, 0,sqr2,0, 1/sqr2,1/sqr2,0};
    memcpy(coord2,coord4,12*sizeof(double));
    ++nt;
    err += test_sasa(s,ref,tolerance,"Four spheres in plane, rotated 45 deg round z-axis.");

    //rotate 90 degrees round x-axis
    double coord5[12] = {-1./sqr2,0,1/sqr2, 0,0,0, 0,0,sqr2, 1/sqr2,0,1/sqr2};
    memcpy(coord2,coord5,12*sizeof(double));
    ++nt;
    err += test_sasa(s,ref,tolerance,"Four spheres in plane, rotated 90 deg round x-axis.");

    printf("\tPassed %d/%d tests (relative errors < %f)\n",nt-err,nt,tolerance);
    number_tests += nt;
    return err;
}

int main(int argc, char **argv) {
    int n_err = 0;
    number_tests = 0;

    sasalib_t *sr = sasalib_init();
    sasalib_t *lr = sasalib_init();

    sasalib_set_algorithm(sr,SASALIB_SHRAKE_RUPLEY);
    sasalib_set_algorithm(lr,SASALIB_LEE_RICHARDS);

    sasalib_set_sr_points(sr,5000);
    sasalib_set_lr_delta(lr,1e-4);

    n_err += test_sasa_alg(sr,1e-3);
    n_err += test_sasa_alg(lr,1e-5);
    
    printf("Total: there were %d errors.\n",n_err);

    return 0;
}
