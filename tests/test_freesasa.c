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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <check.h>
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include <freesasa.h>

#define PASS 1
#define NOPASS 0

double total_ref, polar_ref, apolar_ref;
double tolerance;

freesasa_parameters parameters;
freesasa_result result;

extern int freesasa_fail(const char*);
extern int freesasa_warn(const char*);

double rel_err(double v1, double v2) {
    return fabs(v1-v2)/(fabs(v1)+fabs(v2));
}

double surface_hidden_sphere_intersection(double r1, double r2, double d)
{
    if (d > r1 + r2) return 0;
    if (r1+d < r2) return 4*M_PI*r1*r1;
    if (r2+d < r1) return 4*M_PI*r2*r2;
    return M_PI/d * ( r1*(r2*r2 - (d-r1)*(d-r1))
                    +r2*(r1*r1 - (d-r2)*(d-r2)) ) ;
}
double surface_spheres_intersecting(double r1, double r2, double d)
{
    return 4*M_PI*(r1*r1+r2*r2) - surface_hidden_sphere_intersection(r1,r2,d);
}

double surface_two_spheres(const double *x, const double *r, double probe)
{
    double d2 = (x[0]-x[3])*(x[0]-x[3]) + (x[1]-x[4])*(x[1]-x[4])
        + (x[2]-x[5])*(x[2]-x[5]);
    return surface_spheres_intersecting(r[0]+probe,r[1]+probe,sqrt(d2));
}

int test_sasa(double ref, const char *test, const double *xyz, 
              const double *r, int n)
{
    double err;
    freesasa_result result;
    int pass = PASS;
    freesasa_calc_coord(&result,xyz,r,n,&parameters);
    if ((err = rel_err(ref,result.total))
        > tolerance) {
        pass = NOPASS;
    }
    freesasa_result_free(result);
    return pass;
}

void setup_lr_precision(void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_LEE_RICHARDS;
    parameters.lee_richards_delta = 1e-4;
    tolerance = 1e-5;
}
void teardown_lr_precision(void)
{
    
}
void setup_sr_precision(void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_SHRAKE_RUPLEY;
    parameters.shrake_rupley_n_points = 5000;
    tolerance = 1e-3;
}
void teardown_sr_precision(void)
{
    
}

START_TEST (test_sasa_alg_basic)
{
    // Two spheres, compare with analytic results
    double coord[6] = {0,0,0,2,0,0};
    double r[2] = {1,2};
    double probe = parameters.probe_radius;
    double n = 2;

    ck_assert(test_sasa(surface_two_spheres(coord,r,probe),
                        "Two intersecting spheres along x-axis.",
                        coord,r,n));
    coord[3] = 0; coord[4] = 2;
    ck_assert(test_sasa(surface_two_spheres(coord,r,probe),
                        "Two intersecting spheres along y-axis.",
                        coord,r,n));

    coord[4] = 0; coord[5] = 2;
    ck_assert(test_sasa(surface_two_spheres(coord,r,probe),
                        "Two intersecting spheres along z-axis.",
                        coord,r,n));

    // Four spheres in a plane, all calculations should give similar results
    double coord2[12] = {0,0,0, 1,0,0, 0,1,0, 1,1,0};
    double r2[4] = {1,1,2,1};
    n = 4;
    freesasa_result result;
    freesasa_calc_coord(&result,coord2,r2,4,&parameters);
    double ref = result.total;

    //translate
    for (int i = 0; i < 12; ++i) coord2[i] += 1.;
    ck_assert(test_sasa(ref,"Four spheres in plane, translated",
                        coord2,r2,n));

    //rotate 90 degrees round z-axis
    double coord3[12] = {0,1,0, 0,0,0, 1,1,0, 1,0,0};
    memcpy(coord2,coord3,12*sizeof(double));
    ck_assert(test_sasa(ref,"Four spheres in plane, rotated 90 deg round z-axis.",
                        coord2,r2,n));

    //rotate -45 degrees round z-axis
    double sqr2 = sqrt(2);
    double coord4[12] = {-1./sqr2,1./sqr2,0, 0,0,0, 0,sqr2,0, 1/sqr2,1/sqr2,0};
    memcpy(coord2,coord4,12*sizeof(double));
    ck_assert(test_sasa(ref,"Four spheres in plane, rotated 45 deg round z-axis.",
                        coord2,r2,n));

    //rotate 90 degrees round x-axis
    double coord5[12] = {-1./sqr2,0,1/sqr2, 0,0,0, 0,0,sqr2, 1/sqr2,0,1/sqr2};
    memcpy(coord2,coord5,12*sizeof(double));
    ck_assert(test_sasa(ref,"Four spheres in plane, rotated 90 deg round x-axis.",
                        coord2,r2,n));
    
    freesasa_result_free(result);
}
END_TEST

START_TEST (test_minimal_calc)
{
    freesasa_result result;
    
    double coord[3] = {0,0,0};
    double r[1] = {1.0};

    ck_assert(freesasa_calc_coord(&result,coord,r,1,NULL) == FREESASA_SUCCESS);

    // access areas
    ck_assert(fabs(result.sasa[0] - result.total) < 1e-10);
    ck_assert(fabs(result.total - (4*M_PI*M_PI*2.4*2.4)));
    
    freesasa_result_free(result);
}
END_TEST


void setup_sr (void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_SHRAKE_RUPLEY;
    parameters.shrake_rupley_n_points = 100;
    total_ref = 4759.86096;
    polar_ref = 2232.23039;
    apolar_ref = 2527.63057;
}
void teardown_sr(void)
{

}

void setup_lr (void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_LEE_RICHARDS;
    parameters.lee_richards_delta = 0.25;
    total_ref = 4728.26159;
    polar_ref = 2211.41649;
    apolar_ref = 2516.84510;
}
void teardown_lr(void)
{

}

START_TEST (test_sasa_1ubq)
{
    errno = 0;
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(pdb != NULL);
    freesasa_structure *st = freesasa_structure_from_pdb(pdb,0);
    freesasa_result res;
    double *radii = freesasa_structure_radius(st,NULL);
    ck_assert(freesasa_calc_structure(&res,st,radii,&parameters) == FREESASA_SUCCESS);
    freesasa_strvp* res_class = freesasa_result_classify(res,st,NULL);
    fclose(pdb);

    // The reference values were the output of FreeSASA on 2014-02-10
    ck_assert(fabs(res.total - total_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_POLAR] - polar_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_APOLAR] - apolar_ref) < 1e-5);
    
    FILE *devnull = fopen("/dev/null","w");
    ck_assert(freesasa_write_pdb(devnull,res,st,radii) == FREESASA_SUCCESS);
    ck_assert(freesasa_log(devnull,res,"test",st,NULL,res_class) == FREESASA_SUCCESS);
    ck_assert(freesasa_per_residue_type(devnull,res,st) == FREESASA_SUCCESS);
    ck_assert(freesasa_per_residue(devnull,res,st) == FREESASA_SUCCESS);
    fclose(devnull);

    freesasa_set_verbosity(FREESASA_V_SILENT);
    FILE *nowrite = fopen("/dev/null","r");
    ck_assert(freesasa_log(nowrite,res,"test",st,NULL,res_class) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue_type(nowrite,res,st) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue(nowrite,res,st) == FREESASA_FAIL);
    fclose(nowrite);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    
    free(radii);
    freesasa_strvp_free(res_class);
    freesasa_structure_free(st);
    freesasa_result_free(res);
}
END_TEST

START_TEST (test_trimmed_pdb) 
{
    // This test is due to suggestion from João Rodrigues (issue #6 on Github)
    double total_ref = 15955.547786749;
    double polar_ref = 6543.356616946;
    double apolar_ref = 9412.191169803;
    freesasa_result result;
    freesasa_structure *st;
    FILE *pdb;
    double *radii;
    freesasa_strvp *res_class;
    
    errno = 0;
    pdb = fopen(DATADIR "3bzd_trimmed.pdb","r");
    ck_assert(pdb != NULL);
    st = freesasa_structure_from_pdb(pdb,0);
    fclose(pdb);

    radii = freesasa_structure_radius(st,NULL);
    ck_assert(radii != NULL);
    ck_assert(freesasa_calc_structure(&result,st,radii,NULL) == FREESASA_SUCCESS);
    res_class = freesasa_result_classify(result,st,NULL);
    ck_assert(res_class != NULL);
    
    ck_assert(fabs(result.total - total_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_POLAR] - polar_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_APOLAR] - apolar_ref) < 1e-5);
    
    freesasa_structure_free(st);
    freesasa_result_free(result);
    free(radii);
    freesasa_strvp_free(res_class);
}
END_TEST
/*
START_TEST (test_user_classes)
{
    freesasa *s = freesasa_new(), *s_ref = freesasa_new();
    FILE *config = fopen(DATADIR "test.config", "r");
    ck_assert(config);
    ck_assert(freesasa_user_classification(s,config) == FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_radius(s,"AA","aa") - 1.0) < 1e-5);
    fclose(config);

    config = fopen(DATADIR "oons.config", "r");
    ck_assert(config);
    ck_assert(freesasa_user_classification(s,config) == FREESASA_SUCCESS);
    FILE *ubq = fopen(DATADIR "1ubq.pdb", "r");
    ck_assert(freesasa_calc_pdb(s,ubq) == FREESASA_SUCCESS);
    rewind(ubq);
    freesasa_calc_pdb(s_ref,ubq);    
    ck_assert(fabs(freesasa_area_total(s) - freesasa_area_total(s_ref)) < 1e-5);
    fclose(config);
    


    fclose(ubq);
    freesasa_free(s);
}
END_TEST
*/
/*
START_TEST (test_calc_errors)
{
    freesasa *s = freesasa_new();
    double dummy;

    fputs("Testing error messages:\n",stderr);
    freesasa_fail("Test fail-message.");
    freesasa_warn("Test warn-message.");
    ck_assert(freesasa_get_verbosity() == 0);
    freesasa_set_verbosity(FREESASA_V_SILENT);

    //test empty PDB-file
    FILE *empty = fopen(DATADIR "empty.pdb","r");
    ck_assert(empty != NULL);
    ck_assert(freesasa_calc_pdb(s,empty) == FREESASA_FAIL);
    fclose(empty);

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
    ck_assert(freesasa_radius_atom(s,0) > 0);
    
    freesasa_set_verbosity(0);

    freesasa_free(s);
}
END_TEST

START_TEST (test_multi_calc)
{
#if HAVE_LIBPTHREAD
    errno = 0;
    freesasa *s = freesasa_new();
    //S&R
    freesasa_set_algorithm(s,FREESASA_SHRAKE_RUPLEY);
    freesasa_set_sr_points(s,100);
    freesasa_set_nthreads(s,2);
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
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
*/

Suite *sasa_suite()
{
    Suite *s = suite_create("SASA-calculation");

    TCase *tc_basic = tcase_create("Basic API");
    tcase_add_test(tc_basic, test_minimal_calc);
    //tcase_add_test(tc_basic, test_copyparam);
    //tcase_add_test(tc_basic, test_calc_errors);
    //tcase_add_test(tc_basic, test_user_classes);
    
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

    TCase *tc_trimmed = tcase_create("Trimmed PDB file");
    tcase_add_test(tc_basic, test_trimmed_pdb);

    suite_add_tcase(s, tc_basic);
    suite_add_tcase(s, tc_lr_basic);
    suite_add_tcase(s, tc_sr_basic);
    suite_add_tcase(s, tc_lr);
    suite_add_tcase(s, tc_sr);
    suite_add_tcase(s, tc_trimmed);

#if HAVE_LIBPTHREAD
    printf("Using pthread\n");
    TCase *tc_pthr = tcase_create("Pthread");
    //tcase_add_test(tc_pthr,test_multi_calc);
    suite_add_tcase(s, tc_pthr);
#endif
    return s;
}
