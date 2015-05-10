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

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <check.h>
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <freesasa.h>

#ifndef PI
#define PI 3.14159265358979323846
#endif

#define PASS 1
#define NOPASS 0

freesasa_t *st = NULL;
double total_ref, polar_ref, apolar_ref;
double tolerance;

extern int freesasa_fail(const char*);
extern int freesasa_warn(const char*);
extern int freesasa_set_verbosity(int);
extern int freesasa_get_verbosity();

double rel_err(double v1, double v2) {
    return fabs(v1-v2)/(fabs(v1)+fabs(v2));
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

int test_sasa(double ref, const char *test)
{
    freesasa_refresh(st);
    double err;
    if ((err = rel_err(ref,freesasa_area_total(st)))
        > tolerance) {
        return NOPASS;
    }
    return PASS;
}

void setup_lr_precision(void)
{
    if (st) freesasa_free(st);
    st = freesasa_init();
    freesasa_set_algorithm(st,FREESASA_LEE_RICHARDS);
    freesasa_set_lr_delta(st,1e-4);
    tolerance = 1e-5;
}
void teardown_lr_precision(void)
{
    freesasa_free(st);
}
void setup_sr_precision(void)
{
    if (st) freesasa_free(st);
    st = freesasa_init();
    freesasa_set_algorithm(st,FREESASA_SHRAKE_RUPLEY);
    freesasa_set_sr_points(st,5000);
    tolerance = 1e-3;
}
void teardown_sr_precision(void)
{
    freesasa_free(st);
}

START_TEST (test_sasa_alg_basic)
{
    // Two spheres, compare with analytic results
    double coord[6] = {0,0,0,2,0,0};
    double r[2] = {1,2};
    double probe = freesasa_get_probe_radius(st);
    freesasa_link_coord(st,coord,r,2);

    ck_assert(test_sasa(surface_two_spheres(coord,r,probe),
                        "Two intersecting spheres along x-axis."));
    coord[3] = 0; coord[4] = 2;
    ck_assert(test_sasa(surface_two_spheres(coord,r,probe),
                        "Two intersecting spheres along y-axis."));
    coord[4] = 0; coord[5] = 2;

    ck_assert(test_sasa(surface_two_spheres(coord,r,probe),
                        "Two intersecting spheres along z-axis."));

    // Four spheres in a plane, all calculations should give similar results
    double coord2[12] = {0,0,0, 1,0,0, 0,1,0, 1,1,0};
    double r2[4] = {1,1,2,1};
    freesasa_link_coord(st,coord2,r2,4);
    freesasa_refresh(st);
    double ref = freesasa_area_total(st);

    //translate
    for (int i = 0; i < 12; ++i) coord2[i] += 1.;
    ck_assert(test_sasa(ref,"Four spheres in plane, translated"));

    //rotate 90 degrees round z-axis
    double coord3[12] = {0,1,0, 0,0,0, 1,1,0, 1,0,0};
    memcpy(coord2,coord3,12*sizeof(double));
    ck_assert(test_sasa(ref,"Four spheres in plane, rotated 90 deg round z-axis."));

    //rotate -45 degrees round z-axis
    double sqr2 = sqrt(2);
    double coord4[12] = {-1./sqr2,1./sqr2,0, 0,0,0, 0,sqr2,0, 1/sqr2,1/sqr2,0};
    memcpy(coord2,coord4,12*sizeof(double));
    ck_assert(test_sasa(ref,"Four spheres in plane, rotated 45 deg round z-axis."));

    //rotate 90 degrees round x-axis
    double coord5[12] = {-1./sqr2,0,1/sqr2, 0,0,0, 0,0,sqr2, 1/sqr2,0,1/sqr2};
    memcpy(coord2,coord5,12*sizeof(double));
    ck_assert(test_sasa(ref,"Four spheres in plane, rotated 90 deg round x-axis."));

}
END_TEST

void setup_sr (void)
{
    if (st) freesasa_free(st);
    st = freesasa_init();
    freesasa_set_algorithm(st,FREESASA_SHRAKE_RUPLEY);
    freesasa_set_sr_points(st,100);
    total_ref = 4759.86096;
    polar_ref = 2232.23039;
    apolar_ref = 2527.63057;
}
void teardown_sr(void)
{
    freesasa_free(st);
}

void setup_lr (void)
{
    if (st) freesasa_free(st);
    st = freesasa_init();
    freesasa_set_algorithm(st,FREESASA_LEE_RICHARDS);
    freesasa_set_lr_delta(st,0.25);
    total_ref = 4728.26159;
    polar_ref = 2211.41649;
    apolar_ref = 2516.84510;
}
void teardown_lr(void)
{
    freesasa_free(st);
}

START_TEST (test_sasa_1ubq)
{
    errno = 0;
    FILE *pdb = fopen("data/1ubq.pdb","r");
    if (pdb == NULL) {
        fprintf(stderr,"error reading PDB-file for test. "
                "(Tests must be run from directory test/): %s\n",
                strerror(errno));
    }
    ck_assert(pdb != NULL);
    ck_assert(freesasa_calc_pdb(st,pdb) == FREESASA_SUCCESS);
    fclose(pdb);

    // The reference values were the output of FreeSASA on 2014-02-10
    ck_assert(fabs(freesasa_area_total(st) - total_ref) < 1e-5);
    ck_assert(fabs(freesasa_area_class(st,FREESASA_POLAR) - polar_ref) < 1e-5);
    ck_assert(fabs(freesasa_area_class(st,FREESASA_APOLAR) - apolar_ref) < 1e-5);
    ck_assert(freesasa_area_residue(st,"ALA") > 0);

    const double *r = freesasa_radius_atom_array(st);
    ck_assert(r != NULL);
    ck_assert(r[0] > 0);
    ck_assert(fabs(r[0] - freesasa_radius_atom(st,0)) < 1e-5);

    FILE *devnull = fopen("/dev/null","w");
    ck_assert(freesasa_write_pdb(devnull,st) == FREESASA_SUCCESS);
    ck_assert(freesasa_log(devnull,st) == FREESASA_SUCCESS);
    ck_assert(freesasa_per_residue_type(devnull,st) == FREESASA_SUCCESS);
    fclose(devnull);
}
END_TEST

START_TEST (test_freesasa_api_basic)
{
    extern int freesasa_set_verbosity();
    freesasa_set_verbosity(1);
    freesasa_t *s = freesasa_init();
    ck_assert(s != NULL);
    ck_assert(freesasa_n_atoms(s) == 0);

    // algorithm
    ck_assert(freesasa_algorithm_name(s) != NULL);
    ck_assert(freesasa_set_algorithm(s,-1) == FREESASA_WARN);
    ck_assert(freesasa_set_algorithm(s,1000) == FREESASA_WARN);
    ck_assert(freesasa_set_algorithm(s,FREESASA_LEE_RICHARDS) == FREESASA_SUCCESS);
    ck_assert(freesasa_get_algorithm(s) == FREESASA_LEE_RICHARDS);
    ck_assert(freesasa_algorithm_name(s) != NULL);

    // atom radius calculation (more extensive analysis in test_classify)
    ck_assert(fabs(freesasa_radius("ALA"," H  ")) < 1e-10);

    // probe_radius
    ck_assert(freesasa_set_probe_radius(s,-1.) == FREESASA_WARN);
    ck_assert(freesasa_set_probe_radius(s,1.2) == FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_get_probe_radius(s)-1.2) < 1e-10);

    // L&R delta
    double lrd_def = freesasa_get_lr_delta(s);
    ck_assert(freesasa_set_lr_delta(s,0.5) == FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_get_lr_delta(s)-0.5) < 1e-10);
    ck_assert(freesasa_set_lr_delta(s, 10) == FREESASA_WARN);
    ck_assert(fabs(freesasa_get_lr_delta(s)-10) < 1e-10);
    ck_assert(freesasa_set_lr_delta(s,-1.0) == FREESASA_WARN);
    ck_assert(fabs(freesasa_get_lr_delta(s)-lrd_def) < 1e-10);
    ck_assert(freesasa_get_sr_points(s) == FREESASA_WARN);

    // S&R test-points
    ck_assert(freesasa_set_algorithm(s,FREESASA_SHRAKE_RUPLEY) == FREESASA_SUCCESS);
    int srp_def = freesasa_get_sr_points(s);
    ck_assert(freesasa_set_sr_points(s,100) == FREESASA_SUCCESS);
    ck_assert(freesasa_get_sr_points(s) == 100);
    ck_assert(freesasa_set_sr_points(s,1123) == FREESASA_WARN);
    ck_assert(freesasa_set_sr_points(s,-1123) == FREESASA_WARN);
    ck_assert(freesasa_get_sr_points(s) == srp_def);
    ck_assert(freesasa_get_lr_delta(s) < 0);

    // names
    freesasa_set_proteinname(s,"bla");
    ck_assert_str_eq(freesasa_get_proteinname(s),"bla");

#if HAVE_LIBPTHREAD
    // Threads
    int nt_def = freesasa_get_nthreads(s);
    ck_assert(freesasa_set_nthreads(s,2) == FREESASA_SUCCESS);
    ck_assert(freesasa_get_nthreads(s) == 2);
    ck_assert(freesasa_set_nthreads(s,-1) == FREESASA_WARN);
    ck_assert(freesasa_get_nthreads(s) == nt_def);
#endif

    // Check that results cannot be accessed before calculations are
    // performed
    ck_assert(freesasa_area_total(s) < 0);
    ck_assert(freesasa_area_class(s, FREESASA_POLAR) < 0);
    ck_assert(freesasa_area_class(s, FREESASA_APOLAR) < 0);
    ck_assert(freesasa_per_residue(stdout,s) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue(NULL,s) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue_type(stdout,s) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue_type(NULL,s) == FREESASA_FAIL);
    ck_assert(freesasa_area_residue(s,"ALA") < 0);
    ck_assert(freesasa_write_pdb(stdout,s) == FREESASA_FAIL);
    ck_assert(freesasa_area_atom(s,0) < 0);
    ck_assert(freesasa_area_atom_array(s) == NULL);

    ck_assert(freesasa_log(stdout,s) == FREESASA_WARN);

    freesasa_set_verbosity(0);
}
END_TEST

START_TEST (test_copyparam)
{
    int alg = FREESASA_LEE_RICHARDS;
    int n_sr = 500;
    int n_thread = 3;
    double d_lr = 0.5;
    double probe_radius = 2;
    freesasa_t *s1 = freesasa_init(), *s2 = freesasa_init();
    ck_assert(s1 != NULL);
    ck_assert(s2 != NULL);
    freesasa_set_sr_points(s1,n_sr);
    freesasa_set_algorithm(s1,alg);
    freesasa_set_lr_delta(s1,d_lr);
    freesasa_set_nthreads(s1,n_thread);
    freesasa_set_probe_radius(s1,probe_radius);
    freesasa_copy_param(s2,s1);
    ck_assert(freesasa_get_algorithm(s2) == alg);
    ck_assert(freesasa_get_sr_points(s2) == FREESASA_WARN);
    ck_assert(freesasa_get_nthreads(s2) == n_thread);
    ck_assert(fabs(freesasa_get_lr_delta(s2) - d_lr) < 1e-10);
    ck_assert(fabs(freesasa_get_probe_radius(s2) - probe_radius) < 1e-10);
    freesasa_set_algorithm(s2,FREESASA_SHRAKE_RUPLEY);
    ck_assert(freesasa_get_sr_points(s2) == n_sr);
    freesasa_free(s1);
    freesasa_free(s2);
}
END_TEST

START_TEST (test_minimal_calc)
{
    freesasa_t *s = freesasa_init();

    double coord[3] = {0,0,0};
    double r[1] = {1.0};

    freesasa_set_verbosity(1);

    ck_assert(freesasa_calc_coord(s,coord,r,1) == FREESASA_SUCCESS);

    // access areas
    ck_assert(freesasa_area_atom(s,-1) < 0);
    ck_assert(freesasa_area_atom(s,1) < 0);
    ck_assert(fabs(freesasa_area_atom(s,0) - freesasa_area_total(s)) < 1e-10);
    const double *a = freesasa_area_atom_array(s);
    ck_assert(a != NULL);
    ck_assert(fabs(a[0] - freesasa_area_total(s)) < 1e-10);

    // radii should not have been stored, verify
    ck_assert(freesasa_radius_atom(s,1) < 0);
    ck_assert(freesasa_radius_atom(s,-1) < 0);
    ck_assert(freesasa_radius_atom(s,0) < 0);
    ck_assert(freesasa_radius_atom_array(s) == NULL);

    freesasa_free(s);
    freesasa_set_verbosity(0);
}
END_TEST

START_TEST (test_calc_errors)
{
    freesasa_t *s = freesasa_init();
    double dummy;

    fputs("Testing error messages:\n",stderr);
    freesasa_fail("Test fail-message.");
    freesasa_warn("Test warn-message.");
    ck_assert(freesasa_get_verbosity() == 0);
    freesasa_set_verbosity(1);

    //test empty coordinates
    ck_assert(freesasa_calc_coord(s,&dummy,NULL,0) == FREESASA_WARN);
    ck_assert(freesasa_radius_atom(s,0) == FREESASA_FAIL);
    ck_assert(freesasa_radius_atom_array(s) == NULL);
    ck_assert(freesasa_calc_coord(s,&dummy,&dummy,0) == FREESASA_WARN);
    freesasa_set_algorithm(s,FREESASA_LEE_RICHARDS);
    ck_assert(freesasa_calc_coord(s,&dummy,&dummy,0) == FREESASA_WARN);
    ck_assert(freesasa_radius_atom(s,0) == FREESASA_FAIL);

    //test empty PDB-file
    FILE *empty = fopen("data/empty.pdb","r");
    ck_assert(empty != NULL);
    ck_assert(freesasa_calc_pdb(s,empty) == FREESASA_FAIL);
    fclose(empty);

    //test refresh
    ck_assert(freesasa_refresh(s) == FREESASA_FAIL);
    ck_assert(freesasa_link_coord(s,&dummy,NULL,0) == FREESASA_SUCCESS);
    ck_assert(freesasa_refresh(s) == FREESASA_FAIL);
    freesasa_link_coord(s,&dummy,&dummy,0);
    ck_assert(freesasa_refresh(s) == FREESASA_WARN);

    //test freesasa_calc_atoms
    const char *da[1] = {""};
    const char *ala[1] = {"ALA"};
    const char *calpha[1] = {" CA "};
    const char *unk[1] = {" XY "};
    const char *ill[1] = {"  XY  "};
    const double coord[3] = {0,0,0};

    ck_assert(freesasa_calc_atoms(s,coord,da,da,1) == FREESASA_WARN);
    ck_assert(freesasa_calc_atoms(s,coord,ala,da,1) == FREESASA_WARN);
    ck_assert(freesasa_calc_atoms(s,coord,da,calpha,1) == FREESASA_WARN);
    ck_assert(freesasa_calc_atoms(s,coord,ala,unk,1) == FREESASA_WARN);
    ck_assert(freesasa_calc_atoms(s,coord,ala,ill,1) == FREESASA_WARN);
    ck_assert(freesasa_calc_atoms(s,coord,ill,calpha,1) == FREESASA_WARN);
    ck_assert(freesasa_calc_atoms(s,coord,ala,calpha,1) == FREESASA_SUCCESS);
    ck_assert((dummy=freesasa_area_total(s)) > 0);
    ck_assert(fabs(freesasa_area_class(s,FREESASA_APOLAR)- dummy) < 1e-10);
    ck_assert(fabs(freesasa_area_class(s,FREESASA_POLAR)) < 1e-10 );

    freesasa_set_verbosity(0);

    freesasa_free(s);
}
END_TEST

START_TEST (test_multi_calc)
{
#if HAVE_LIBPTHREAD
    errno = 0;
    freesasa_t *s = freesasa_init();
    //S&R
    freesasa_set_algorithm(s,FREESASA_SHRAKE_RUPLEY);
    freesasa_set_sr_points(s,100);
    freesasa_set_nthreads(s,2);
    FILE *pdb = fopen("data/1ubq.pdb","r");
    if (pdb == NULL) {
        fprintf(stderr,"error reading PDB-file for test. "
                "(Tests must be run from directory test/): %s\n",
                strerror(errno));
    }
    ck_assert(pdb != NULL);
    ck_assert(freesasa_calc_pdb(s,pdb) == FREESASA_SUCCESS);
    // The reference values were the output of FreeSASA on 2014-02-10
    ck_assert(fabs(freesasa_area_total(s) - 4759.86096) < 1e-5);
    // L&R
    freesasa_set_algorithm(s,FREESASA_LEE_RICHARDS);
    freesasa_set_lr_delta(s,0.25);
    rewind(pdb);
    ck_assert(freesasa_calc_pdb(s,pdb) == FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_area_total(s) - 4728.26159) < 1e-5);

    fclose(pdb);
    freesasa_free(s);
#endif
}
END_TEST


Suite *sasa_suite()
{
    Suite *s = suite_create("SASA-calculation");

    TCase *tc_basic = tcase_create("Basic API");
    tcase_add_test(tc_basic, test_freesasa_api_basic);
    tcase_add_test(tc_basic, test_copyparam);
    tcase_add_test(tc_basic, test_calc_errors);
    tcase_add_test(tc_basic, test_minimal_calc);

    TCase *tc_lr_basic = tcase_create("Basic L&R");
    tcase_add_checked_fixture(tc_lr_basic,setup_lr_precision,teardown_lr_precision);
    tcase_add_test(tc_lr_basic, test_sasa_alg_basic);

    TCase *tc_sr_basic = tcase_create("Basic S&R");
    tcase_add_checked_fixture(tc_sr_basic,setup_sr_precision,teardown_sr_precision);
    tcase_add_test(tc_sr_basic, test_sasa_alg_basic);

    TCase *tc_lr = tcase_create("1UBQ-L&R");
    tcase_add_checked_fixture(tc_lr,setup_lr,teardown_lr);
    tcase_add_test(tc_lr, test_sasa_1ubq);

    TCase *tc_sr = tcase_create("1UBQ-S&R");
    tcase_add_checked_fixture(tc_sr,setup_sr,teardown_sr);
    tcase_add_test(tc_sr, test_sasa_1ubq);

    suite_add_tcase(s, tc_basic);
    suite_add_tcase(s, tc_lr_basic);
    suite_add_tcase(s, tc_sr_basic);
    suite_add_tcase(s, tc_lr);
    suite_add_tcase(s, tc_sr);

#if HAVE_LIBPTHREAD
    printf("Using pthread\n");
    TCase *tc_pthr = tcase_create("Pthread");
    tcase_add_test(tc_pthr,test_multi_calc);
    suite_add_tcase(s, tc_pthr);
#endif
    return s;
}
