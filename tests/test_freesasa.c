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
#include <structure.h>
#include "tools.h"

#define PASS 1
#define NOPASS 0

double total_ref, polar_ref, apolar_ref;
double tolerance;

freesasa_parameters parameters;

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
    freesasa_result *result;
    int pass = PASS;
    result = freesasa_calc_coord(xyz,r,n,&parameters);
    if ((err = rel_err(ref,result->total))
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
    parameters.lee_richards_n_slices = 20000;
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
    freesasa_result *result = freesasa_calc_coord(coord2,r2,4,&parameters);
    double ref = result->total;
    freesasa_result_free(result);

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
    
}
END_TEST

START_TEST (test_minimal_calc)
{
    freesasa_result *result;
    
    double coord[3] = {0,0,0};
    double r[1] = {1.0};

    ck_assert((result = freesasa_calc_coord(coord,r,1,NULL)) != NULL);

    // access areas
    ck_assert(fabs(result->sasa[0] - result->total) < 1e-10);
    ck_assert(fabs(result->total - (4*M_PI*M_PI*2.4*2.4)));
    
    freesasa_result_free(result);
}
END_TEST


void setup_sr (void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_SHRAKE_RUPLEY;
    parameters.shrake_rupley_n_points = 100;
    parameters.n_threads = 1;
    total_ref = 4779.5109924;
    polar_ref = 2236.9298941;
    apolar_ref = 2542.5810983;
}
void teardown_sr(void)
{

}

void setup_lr (void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_LEE_RICHARDS;
    parameters.lee_richards_n_slices = 20;
    parameters.n_threads = 1;
    total_ref = 4759.46651;
    polar_ref = 2226.83182;
    apolar_ref = 2532.63469;
}
void teardown_lr(void)
{

}

START_TEST (test_sasa_1ubq)
{
    errno = 0;
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(pdb != NULL);
    freesasa_structure *st = freesasa_structure_from_pdb(pdb,NULL,0);
    freesasa_result *res;
    ck_assert((res = freesasa_calc_structure(st,&parameters)) != NULL);
    freesasa_strvp* res_class = freesasa_result_classify(res,st,NULL);
    fclose(pdb);

    // The reference values were the output of FreeSASA on 2014-02-10
    ck_assert(fabs(res->total - total_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_POLAR] - polar_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_APOLAR] - apolar_ref) < 1e-5);
    
    FILE *devnull = fopen("/dev/null","w");
    ck_assert(freesasa_write_pdb(devnull,res,st) == FREESASA_SUCCESS);
    ck_assert(freesasa_log(devnull,res,"test",NULL,res_class) == FREESASA_SUCCESS);
    ck_assert(freesasa_per_residue_type(devnull,res,st) == FREESASA_SUCCESS);
    ck_assert(freesasa_per_residue(devnull,res,st) == FREESASA_SUCCESS);
    fclose(devnull);

    freesasa_set_verbosity(FREESASA_V_SILENT);
    FILE *nowrite = fopen("/dev/null","r");
    ck_assert(freesasa_log(nowrite,res,"test",NULL,res_class) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue_type(nowrite,res,st) == FREESASA_FAIL);
    ck_assert(freesasa_per_residue(nowrite,res,st) == FREESASA_FAIL);
    fclose(nowrite);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    
    freesasa_strvp_free(res_class);
    freesasa_structure_free(st);
    freesasa_result_free(res);
}
END_TEST

START_TEST (test_write_1ubq) {
    FILE *tf = fopen("tmp/dummy_bfactors.pdb","w+"),
        *ref = fopen(DATADIR "reference_bfactors.pdb","r"),
        *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(tf != NULL);
    ck_assert(ref != NULL);
    ck_assert(pdb != NULL);
    freesasa_structure *s = freesasa_structure_from_pdb(pdb, NULL, 0);
    const int n = freesasa_structure_n(s);
    freesasa_result res;
    fclose(pdb);

    res.sasa = malloc(sizeof(double)*n);
    for (int i = 0; i < n; ++i) res.sasa[i] = 1.23;
    
    freesasa_structure_set_radius(s, res.sasa);
    ck_assert(freesasa_write_pdb(tf, &res, s) == FREESASA_SUCCESS);

    rewind(tf);
    free(res.sasa);

    //check that output matches reference file
    size_t bufsize = 100;
    char *buf_tf = malloc(bufsize), *buf_ref = malloc(bufsize);
    while(getline(&buf_tf,&bufsize,tf) > 0 && getline(&buf_ref,&bufsize,ref) > 0) {
        ck_assert_str_eq(buf_ref,buf_tf);
    }
    free(buf_tf);
    free(buf_ref);
    fclose(ref);
    fclose(tf);
}
END_TEST


START_TEST (test_trimmed_pdb) 
{
    // This test is due to suggestion from João Rodrigues (issue #6 on Github)
    double total_ref = 15964.2723037;
    double polar_ref = 6545.4595991;
    double apolar_ref = 9418.8127046;
    freesasa_result *result;
    freesasa_structure *st;
    FILE *pdb;
    freesasa_strvp *res_class;
    
    errno = 0;
    pdb = fopen(DATADIR "3bzd_trimmed.pdb","r");
    ck_assert(pdb != NULL);
    st = freesasa_structure_from_pdb(pdb,NULL,0);
    fclose(pdb);

    ck_assert((result = freesasa_calc_structure(st,NULL)) != NULL);
    res_class = freesasa_result_classify(result,st,NULL);
    ck_assert(res_class != NULL);
    
    ck_assert(fabs(result->total - total_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_POLAR] - polar_ref) < 1e-5);
    ck_assert(fabs(res_class->value[FREESASA_APOLAR] - apolar_ref) < 1e-5);
    
    freesasa_structure_free(st);
    freesasa_result_free(result);
    freesasa_strvp_free(res_class);
}
END_TEST

START_TEST (test_user_classes)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    FILE *config = fopen(DATADIR "oons.config", "r");
    freesasa_structure *st, *st_ref;
    freesasa_classifier *user_classifier;
    freesasa_result *res;
    freesasa_strvp *res_class, *res_class_ref;
    const double *radii, *radii_ref;

    ck_assert(pdb != NULL);
    ck_assert(config != NULL);

    user_classifier = freesasa_classifier_from_file(config);
    ck_assert(user_classifier != NULL);
    fclose(config);

    st = freesasa_structure_from_pdb(pdb, user_classifier, 0);
    ck_assert(st != NULL);
    rewind(pdb);
    st_ref = freesasa_structure_from_pdb(pdb, NULL, 0);
    ck_assert(st_ref != NULL);
    fclose(pdb);
    
    radii = freesasa_structure_radius(st);
    radii_ref = freesasa_structure_radius(st_ref);
    ck_assert(radii != NULL);
    ck_assert(radii_ref != NULL);
    for (int i = 0; i < freesasa_structure_n(st); ++i) {
        ck_assert(fabs(radii[i] - radii_ref[i]) < 1e-5);
    }
    ck_assert((res = freesasa_calc_structure(st,NULL)) != NULL);
    res_class = freesasa_result_classify(res,st,user_classifier);
    res_class_ref = freesasa_result_classify(res,st,NULL);
    ck_assert_int_eq(res_class->n, res_class_ref->n);
    for (int i = 0; i < res_class->n; ++i) {
        ck_assert(fabs(res_class->value[i] - res_class_ref->value[i]) < 1e-10);
    }
    
    freesasa_strvp_free(res_class);
    freesasa_strvp_free(res_class_ref);
    freesasa_structure_free(st);
    freesasa_classifier_free(user_classifier);
    freesasa_result_free(res);
}
END_TEST


START_TEST (test_calc_errors)
{
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    fputs("Testing error messages:\n",stderr);
    freesasa_fail("Test fail-message.");
    freesasa_warn("Test warn-message.");
    ck_assert(freesasa_get_verbosity() == 0);
    freesasa_set_verbosity(FREESASA_V_SILENT);

    //test empty PDB-file
    FILE *empty = fopen(DATADIR "empty.pdb","r");
    ck_assert(empty != NULL);
    ck_assert(freesasa_structure_from_pdb(empty, NULL, 0) == NULL);
    fclose(empty);
    
    freesasa_set_verbosity(FREESASA_V_NORMAL);

}
END_TEST

START_TEST (test_multi_calc)
{
#if HAVE_LIBPTHREAD
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    freesasa_structure *st = freesasa_structure_from_pdb(pdb, NULL, 0);
    freesasa_result *res;
    freesasa_parameters p = freesasa_default_parameters;

    fclose(pdb);

    //S&R
    p.n_threads = 2;
    p.alg = FREESASA_SHRAKE_RUPLEY;
    ck_assert((res = freesasa_calc_structure(st,&p)) != NULL);
    ck_assert(fabs(res->total - 4779.5109924) < 1e-5);
    // L&R
    p.alg = FREESASA_LEE_RICHARDS;
    p.lee_richards_n_slices = 20;
    ck_assert((res = freesasa_calc_structure(st,&p)) != NULL);
    ck_assert(fabs(res->total - 4759.46651) < 1e-5);
    
    freesasa_structure_free(st);
    freesasa_result_free(res);
#endif
}
END_TEST

// test an NMR structure with hydrogens and several models
START_TEST (test_1d3z) 
{
    FILE *pdb = fopen(DATADIR "1d3z.pdb","r");
    int n = 0;
    freesasa_structure* st = freesasa_structure_from_pdb(pdb, NULL, 0);
    freesasa_result *result = freesasa_calc_structure(st, NULL);
    double *radii_ref = malloc(sizeof(double)*602);
    ck_assert(freesasa_structure_n(st) == 602);
    ck_assert(fabs(result->total - 4945.8705756) < 1e-5);
    memcpy(radii_ref,freesasa_structure_radius(st),sizeof(double)*602);
    rewind(pdb);
    freesasa_structure_free(st);
    
    freesasa_set_verbosity(FREESASA_V_SILENT);
    st = freesasa_structure_from_pdb(pdb, NULL, FREESASA_INCLUDE_HYDROGEN);
    result = freesasa_calc_structure(st, NULL);
    ck_assert(freesasa_structure_n(st) == 1231);
    ck_assert(fabs(result->total - 4988.93) < 0.1);
    rewind(pdb);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    freesasa_structure** ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_MODELS);
    ck_assert(n == 10);
    result = freesasa_calc_structure(ss[0], NULL);
    ck_assert(freesasa_structure_n(ss[0]) == 602);
    ck_assert(fabs(result->total - 4945.8705756) < 1e-5);
    for (int i = 0; i < n; ++i) {
        const double *r2 = freesasa_structure_radius(ss[i]);
        ck_assert(r2 != NULL);
        for (int j = 0; j < 602; ++j) {
            ck_assert(fabs(r2[j] - radii_ref[j]) < 1e-10);
        }
        freesasa_structure_free(ss[i]);
    }
    free(ss);
    freesasa_result_free(result);
    freesasa_structure_free(st);
    fclose(pdb);
}
END_TEST

Suite *sasa_suite()
{
    Suite *s = suite_create("SASA-calculation");

    TCase *tc_basic = tcase_create("API");
    tcase_add_test(tc_basic, test_minimal_calc);
    tcase_add_test(tc_basic, test_calc_errors);
    tcase_add_test(tc_basic, test_user_classes);
    tcase_add_test(tc_basic, test_write_1ubq);
    
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
    tcase_add_test(tc_trimmed, test_trimmed_pdb);

    TCase *tc_1d3z = tcase_create("NMR PDB-file 1D3Z (several models, hydrogens)");
    tcase_add_test(tc_1d3z,test_1d3z);

    suite_add_tcase(s, tc_basic);
    suite_add_tcase(s, tc_lr_basic);
    suite_add_tcase(s, tc_sr_basic);
    suite_add_tcase(s, tc_lr);
    suite_add_tcase(s, tc_sr);
    suite_add_tcase(s, tc_trimmed);
    suite_add_tcase(s, tc_1d3z);

#if HAVE_LIBPTHREAD
    printf("Using pthread\n");
    TCase *tc_pthr = tcase_create("Pthread");
    tcase_add_test(tc_pthr,test_multi_calc);
    suite_add_tcase(s, tc_pthr);
#endif
    return s;
}
