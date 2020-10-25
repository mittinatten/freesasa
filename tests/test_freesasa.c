#include <check.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <freesasa.h>
#include <freesasa_internal.h>

#include "tools.h"

#define PASS 1
#define NOPASS 0

double total_ref, polar_ref, apolar_ref;
double tolerance;

freesasa_parameters parameters;

double rel_err(double v1, double v2)
{
    return fabs(v1 - v2) / (fabs(v1) + fabs(v2));
}

double surface_hidden_sphere_intersection(double r1, double r2, double d)
{
    if (d > r1 + r2) return 0;
    if (r1 + d < r2) return 4 * M_PI * r1 * r1;
    if (r2 + d < r1) return 4 * M_PI * r2 * r2;
    return M_PI / d * (r1 * (r2 * r2 - (d - r1) * (d - r1)) + r2 * (r1 * r1 - (d - r2) * (d - r2)));
}
double surface_spheres_intersecting(double r1, double r2, double d)
{
    return 4 * M_PI * (r1 * r1 + r2 * r2) - surface_hidden_sphere_intersection(r1, r2, d);
}

double surface_two_spheres(const double *x, const double *r, double probe)
{
    double d2 = (x[0] - x[3]) * (x[0] - x[3]) + (x[1] - x[4]) * (x[1] - x[4]) + (x[2] - x[5]) * (x[2] - x[5]);
    return surface_spheres_intersecting(r[0] + probe, r[1] + probe, sqrt(d2));
}

int test_sasa(double ref, const char *test, const double *xyz,
              const double *r, int n)
{
    double err;
    freesasa_result *result;
    int pass = PASS;
    result = freesasa_calc_coord(xyz, r, n, &parameters);
    if ((err = rel_err(ref, result->total)) > tolerance) {
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

START_TEST(test_sasa_alg_basic)
{
    // Two spheres, compare with analytic results
    double coord[6] = {0, 0, 0, 2, 0, 0};
    double r[2] = {1, 2};
    double probe = parameters.probe_radius;
    double n = 2;

    ck_assert(test_sasa(surface_two_spheres(coord, r, probe),
                        "Two intersecting spheres along x-axis.",
                        coord, r, n));
    coord[3] = 0;
    coord[4] = 2;
    ck_assert(test_sasa(surface_two_spheres(coord, r, probe),
                        "Two intersecting spheres along y-axis.",
                        coord, r, n));

    coord[4] = 0;
    coord[5] = 2;
    ck_assert(test_sasa(surface_two_spheres(coord, r, probe),
                        "Two intersecting spheres along z-axis.",
                        coord, r, n));

    // Four spheres in a plane, all calculations should give similar results
    double coord2[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0};
    double r2[4] = {1, 1, 2, 1};
    n = 4;
    freesasa_result *result = freesasa_calc_coord(coord2, r2, 4, &parameters);
    double ref = result->total;
    freesasa_result_free(result);

    //translate
    for (int i = 0; i < 12; ++i)
        coord2[i] += 1.;
    ck_assert(test_sasa(ref, "Four spheres in plane, translated",
                        coord2, r2, n));

    //rotate 90 degrees round z-axis
    double coord3[12] = {0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0};
    memcpy(coord2, coord3, 12 * sizeof(double));
    ck_assert(test_sasa(ref, "Four spheres in plane, rotated 90 deg round z-axis.",
                        coord2, r2, n));

    //rotate -45 degrees round z-axis
    double sqr2 = sqrt(2);
    double coord4[12] = {-1. / sqr2, 1. / sqr2, 0, 0, 0, 0, 0, sqr2, 0, 1 / sqr2, 1 / sqr2, 0};
    memcpy(coord2, coord4, 12 * sizeof(double));
    ck_assert(test_sasa(ref, "Four spheres in plane, rotated 45 deg round z-axis.",
                        coord2, r2, n));

    //rotate 90 degrees round x-axis
    double coord5[12] = {-1. / sqr2, 0, 1 / sqr2, 0, 0, 0, 0, 0, sqr2, 1 / sqr2, 0, 1 / sqr2};
    memcpy(coord2, coord5, 12 * sizeof(double));
    ck_assert(test_sasa(ref, "Four spheres in plane, rotated 90 deg round x-axis.",
                        coord2, r2, n));
}
END_TEST

START_TEST(test_minimal_calc)
{
    double coord[3] = {0, 0, 0};
    double r[1] = {1.0};

    freesasa_result *result = freesasa_calc_coord(coord, r, 1, NULL);

    ck_assert_ptr_ne(result, NULL);

    // access areas
    ck_assert(fabs(result->sasa[0] - result->total) < 1e-10);
    ck_assert(fabs(result->total - (4 * M_PI * M_PI * 2.4 * 2.4)));

    freesasa_result_free(result);
}
END_TEST

void setup_sr(void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_SHRAKE_RUPLEY;
    parameters.shrake_rupley_n_points = 100;
    parameters.n_threads = 1;
    total_ref = 4834.716265;
    polar_ref = 2515.821238;
    apolar_ref = 2318.895027;
}
void teardown_sr(void)
{
}

void setup_lr(void)
{
    parameters = freesasa_default_parameters;
    parameters.alg = FREESASA_LEE_RICHARDS;
    parameters.lee_richards_n_slices = 20;
    parameters.n_threads = 1;
    total_ref = 4804.055641;
    polar_ref = 2504.217302;
    apolar_ref = 2299.838339;
}
void teardown_lr(void)
{
}

START_TEST(test_sasa_1ubq)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    ck_assert(pdb != NULL);
    const freesasa_classifier *classifier = &freesasa_default_classifier;
    freesasa_structure *st = freesasa_structure_from_pdb(pdb, classifier, 0);
    freesasa_result *res;
    freesasa_node *tree = freesasa_tree_new();
    ck_assert((res = freesasa_calc_structure(st, &parameters)) != NULL);

    freesasa_nodearea res_class = freesasa_result_classes(st, res);

    freesasa_tree_add_result(tree, res, st, "test");
    fclose(pdb);

    ck_assert(float_eq(res->total, total_ref, 1e-5));
    ck_assert(float_eq(res_class.polar, polar_ref, 1e-5));
    ck_assert(float_eq(res_class.apolar, apolar_ref, 1e-5));

    FILE *devnull = fopen("/dev/null", "w");
    ck_assert(freesasa_write_pdb(devnull, tree) == FREESASA_SUCCESS);
    ck_assert(freesasa_tree_export(devnull, tree, FREESASA_PDB) == FREESASA_SUCCESS);
    ck_assert(freesasa_tree_export(devnull, tree, FREESASA_LOG) == FREESASA_SUCCESS);
    ck_assert(freesasa_tree_export(devnull, tree, FREESASA_RES) == FREESASA_SUCCESS);
    ck_assert(freesasa_tree_export(devnull, tree, FREESASA_SEQ) == FREESASA_SUCCESS);
    ck_assert(freesasa_tree_export(devnull, tree, FREESASA_RSA) == FREESASA_SUCCESS);
    if (USE_JSON) {
        ck_assert(freesasa_tree_export(devnull, tree, FREESASA_JSON) == FREESASA_SUCCESS);
    } else {
        ck_assert(freesasa_tree_export(devnull, tree, FREESASA_JSON) == FREESASA_FAIL);
    }
    if (USE_XML) {
        ck_assert(freesasa_tree_export(devnull, tree, FREESASA_XML) == FREESASA_SUCCESS);
    } else {
        ck_assert(freesasa_tree_export(devnull, tree, FREESASA_XML) == FREESASA_FAIL);
    }
    fclose(devnull);

    freesasa_set_verbosity(FREESASA_V_SILENT);
    FILE *nowrite = fopen("/dev/null", "r");
    ck_assert(freesasa_write_log(nowrite, tree) == FREESASA_FAIL);
    ck_assert(freesasa_tree_export(nowrite, tree, FREESASA_RSA) == FREESASA_FAIL);
    ck_assert(freesasa_tree_export(nowrite, tree, FREESASA_JSON) == FREESASA_FAIL);
    if (USE_JSON) {
        ck_assert(freesasa_tree_export(nowrite, tree, FREESASA_JSON) == FREESASA_FAIL);
    }
    if (USE_XML) {
        ck_assert(freesasa_tree_export(nowrite, tree, FREESASA_XML) == FREESASA_FAIL);
    }
    fclose(nowrite);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    freesasa_structure_free(st);
    freesasa_result_free(res);
    freesasa_node_free(tree);
}
END_TEST

START_TEST(test_write_pdb)
{
    FILE *tf = fopen("tmp/dummy_bfactors.pdb", "w+"),
         *ref = fopen(DATADIR "reference_bfactors.pdb", "r"),
         *pdb = fopen(DATADIR "1ubq.pdb", "r"),
         *devnull = fopen("/dev/null", "w");
    ck_assert(tf != NULL);
    ck_assert(ref != NULL);
    ck_assert(pdb != NULL);
    freesasa_structure *s = freesasa_structure_from_pdb(pdb, NULL, 0);
    const int n = freesasa_structure_n(s);
    freesasa_result res;
    freesasa_node *root;
    fclose(pdb);

    res.sasa = malloc(sizeof(double) * n);
    for (int i = 0; i < n; ++i)
        res.sasa[i] = 1.23;
    res.parameters = freesasa_default_parameters;
    res.n_atoms = n;

    freesasa_structure_set_radius(s, res.sasa);
    root = freesasa_tree_init(&res, s, "bla");
    ck_assert(freesasa_write_pdb(tf, root) == FREESASA_SUCCESS);

    rewind(tf);
    free(res.sasa);

    //check that output matches reference file
    size_t bufsize = 100;
    char *buf_tf = malloc(bufsize), *buf_ref = malloc(bufsize);
    while (getline(&buf_tf, &bufsize, tf) > 0 && getline(&buf_ref, &bufsize, ref) > 0) {
        // skip remarks
        while (strncmp(buf_ref, "REMARK", 6) == 0 && getline(&buf_ref, &bufsize, ref) > 0)
            ;
        while (strncmp(buf_tf, "REMARK", 6) == 0 && getline(&buf_tf, &bufsize, tf) > 0)
            ;
        ck_assert_str_eq(buf_ref, buf_tf);
    }
    free(buf_tf);
    free(buf_ref);

    freesasa_structure_free(s);
    freesasa_node_free(root);

    // Can't write pdb from structure not initialized from pdb
    s = freesasa_structure_new();
    freesasa_structure_add_atom(s, "C", "ALA", "   1", 'C', 0, 0, 0);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    root = freesasa_tree_init(&res, s, "bla");
    ck_assert_int_eq(freesasa_write_pdb(devnull, root), FREESASA_FAIL);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    freesasa_structure_free(s);
    freesasa_node_free(root);
    fclose(devnull);
    fclose(ref);
    fclose(tf);
}
END_TEST

START_TEST(test_trimmed_pdb)
{
    // This test is due to suggestion from João Rodrigues (issue #6 on Github)
    double total_ref = 16133.867124;
    double polar_ref = 7432.608118;
    double apolar_ref = 8701.259006;
    freesasa_parameters param = freesasa_default_parameters;
    const freesasa_classifier *classifier = &freesasa_default_classifier;
    freesasa_result *result;
    freesasa_structure *st;
    freesasa_nodearea res_class;
    FILE *pdb;
    param.alg = FREESASA_SHRAKE_RUPLEY;

    pdb = fopen(DATADIR "3bzd_trimmed.pdb", "r");
    ck_assert(pdb != NULL);
    st = freesasa_structure_from_pdb(pdb, classifier, 0);
    fclose(pdb);

    result = freesasa_calc_structure(st, &param);
    ck_assert_ptr_ne(result, NULL);
    res_class = freesasa_result_classes(st, result);

    ck_assert(float_eq(result->total, total_ref, 1e-5));
    ck_assert(float_eq(res_class.polar, polar_ref, 1e-5));
    ck_assert(float_eq(res_class.apolar, apolar_ref, 1e-5));

    freesasa_structure_free(st);
    freesasa_result_free(result);
}
END_TEST

START_TEST(test_user_classes)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r"),
         *clf = fopen(SHAREDIR "protor.config", "r");
    freesasa_structure *st, *st_ref;
    freesasa_classifier *user_classifier;
    freesasa_result *res, *res_ref;
    freesasa_nodearea res_class, res_class_ref;
    const double *radii, *radii_ref;

    ck_assert(pdb != NULL);
    ck_assert(clf != NULL);

    user_classifier = freesasa_classifier_from_file(clf);
    fclose(clf);
    ck_assert(user_classifier != NULL);

    st = freesasa_structure_from_pdb(pdb, user_classifier, 0);
    ck_assert(st != NULL);
    rewind(pdb);
    st_ref = freesasa_structure_from_pdb(pdb, &freesasa_protor_classifier, 0);
    ck_assert(st_ref != NULL);
    fclose(pdb);

    radii = freesasa_structure_radius(st);
    radii_ref = freesasa_structure_radius(st_ref);
    ck_assert(radii != NULL);
    ck_assert(radii_ref != NULL);
    for (int i = 0; i < freesasa_structure_n(st); ++i) {
        ck_assert(float_eq(radii[i], radii_ref[i], 1e-10));
    }
    res = freesasa_calc_structure(st, NULL);
    res_ref = freesasa_calc_structure(st_ref, NULL);
    ck_assert_ptr_ne(res, NULL);
    ck_assert_ptr_ne(res, NULL);
    res_class = freesasa_result_classes(st, res);
    res_class_ref = freesasa_result_classes(st_ref, res_ref);
    ck_assert(float_eq(res_class.total, res_class_ref.total, 1e-10));
    ck_assert(float_eq(res_class.polar, res_class_ref.polar, 1e-10));
    ck_assert(float_eq(res_class.apolar, res_class_ref.apolar, 1e-10));
    ck_assert(float_eq(res_class.unknown, res_class_ref.unknown, 1e-10));
    ck_assert(float_eq(res_class.side_chain, res_class_ref.side_chain, 1e-10));
    ck_assert(float_eq(res_class.main_chain, res_class_ref.main_chain, 1e-10));

    freesasa_structure_free(st);
    freesasa_classifier_free(user_classifier);
    freesasa_result_free(res);
}
END_TEST

START_TEST(test_calc_errors)
{
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    fputs("Testing error messages:\n", stderr);
    freesasa_fail("Test fail-message.");
    freesasa_warn("Test warn-message.");
    ck_assert(freesasa_get_verbosity() == 0);
    freesasa_set_verbosity(FREESASA_V_SILENT);

    //test empty PDB-file
    FILE *empty = fopen(DATADIR "empty.pdb", "r");
    ck_assert(empty != NULL);
    ck_assert(freesasa_structure_from_pdb(empty, NULL, 0) == NULL);
    fclose(empty);

    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST(test_multi_calc)
{
#if USE_THREADS
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *st = freesasa_structure_from_pdb(pdb, NULL, 0);
    freesasa_result *res;
    freesasa_parameters p = freesasa_default_parameters;

    fclose(pdb);

    //S&R
    p.n_threads = 2;
    p.alg = FREESASA_SHRAKE_RUPLEY;
    ck_assert((res = freesasa_calc_structure(st, &p)) != NULL);
    ck_assert(fabs(res->total - 4834.716265) < 1e-5);
    // L&R
    p.alg = FREESASA_LEE_RICHARDS;
    p.lee_richards_n_slices = 20;
    ck_assert((res = freesasa_calc_structure(st, &p)) != NULL);
    ck_assert(fabs(res->total - 4804.055641) < 1e-5);

    freesasa_structure_free(st);
    freesasa_result_free(res);
#endif /* USE_THREADS */
}
END_TEST

// test an NMR structure with hydrogens and several models
START_TEST(test_1d3z)
{
    FILE *pdb = fopen(DATADIR "1d3z.pdb", "r");
    int n = 0;
    freesasa_parameters param = freesasa_default_parameters;
    param.alg = FREESASA_SHRAKE_RUPLEY;
    freesasa_structure *st = freesasa_structure_from_pdb(pdb, NULL, 0);
    freesasa_result *result = freesasa_calc_structure(st, &param);
    double *radii_ref = malloc(sizeof(double) * 602);
    ck_assert(freesasa_structure_n(st) == 602);
    ck_assert(fabs(result->total - 5000.340175) < 1e-5);
    memcpy(radii_ref, freesasa_structure_radius(st), sizeof(double) * 602);
    rewind(pdb);
    freesasa_structure_free(st);

    freesasa_set_verbosity(FREESASA_V_SILENT);
    st = freesasa_structure_from_pdb(pdb, NULL, FREESASA_INCLUDE_HYDROGEN);
    result = freesasa_calc_structure(st, &param);
    ck_assert(freesasa_structure_n(st) == 1231);
    ck_assert(fabs(result->total - 5035.614493) < 1e-5);
    rewind(pdb);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    freesasa_structure **ss = freesasa_structure_array(pdb, &n, NULL, FREESASA_SEPARATE_MODELS);
    ck_assert(n == 10);
    result = freesasa_calc_structure(ss[0], &param);
    ck_assert(freesasa_structure_n(ss[0]) == 602);
    ck_assert(fabs(result->total - 5000.340175) < 1e-5);
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

START_TEST(test_memerr)
{
    freesasa_parameters p = freesasa_default_parameters;
    double v[18] = {0, 0, 0, 1, 1, 1, -1, 1, -1, 2, 0, -2, 2, 2, 0, -5, 5, 5};
    struct coord_t coord = {.xyz = v, .n = 6, .is_linked = 0};
    const double r[6] = {4, 2, 2, 2, 2, 2};
    void *ptr;
    p.shrake_rupley_n_points = 10; // so the loop below will be fast

    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (int i = 1; i < 45; ++i) {
        p.alg = FREESASA_SHRAKE_RUPLEY;
        set_fail_after(i);
        ptr = freesasa_calc(&coord, r, &p);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
        p.alg = FREESASA_LEE_RICHARDS;
        set_fail_after(i);
        ptr = freesasa_calc(&coord, r, &p);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
    }

    FILE *file = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *s = freesasa_structure_from_pdb(file, NULL, 0);
    for (int i = 1; i < 256; i *= 2) { //try to spread it out without doing too many calculations
        set_fail_after(i);
        ptr = freesasa_calc_structure(s, NULL);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
        set_fail_after(i);
        ptr = freesasa_structure_get_chains(s, "A", NULL, 0);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
    }
    freesasa_structure_free(s);
    fclose(file);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

extern TCase *test_LR_static();

Suite *sasa_suite()
{
    Suite *s = suite_create("SASA-calculation");

    TCase *tc_basic = tcase_create("API");
    tcase_add_test(tc_basic, test_minimal_calc);
    tcase_add_test(tc_basic, test_calc_errors);
    tcase_add_test(tc_basic, test_user_classes);
    tcase_add_test(tc_basic, test_write_pdb);
    tcase_add_test(tc_basic, test_memerr);

    TCase *tc_lr_basic = tcase_create("Basic L&R");
    tcase_add_checked_fixture(tc_lr_basic, setup_lr_precision, teardown_lr_precision);
    tcase_add_test(tc_lr_basic, test_sasa_alg_basic);

    TCase *tc_lr_static = test_LR_static();

    TCase *tc_sr_basic = tcase_create("Basic S&R");
    tcase_add_checked_fixture(tc_sr_basic, setup_sr_precision, teardown_sr_precision);
    tcase_add_test(tc_sr_basic, test_sasa_alg_basic);

    TCase *tc_lr = tcase_create("1UBQ-L&R");
    tcase_add_checked_fixture(tc_lr, setup_lr, teardown_lr);
    tcase_add_test(tc_lr, test_sasa_1ubq);

    TCase *tc_sr = tcase_create("1UBQ-S&R");
    tcase_add_checked_fixture(tc_sr, setup_sr, teardown_sr);
    tcase_add_test(tc_sr, test_sasa_1ubq);

    TCase *tc_trimmed = tcase_create("Trimmed PDB file");
    tcase_add_test(tc_trimmed, test_trimmed_pdb);

    TCase *tc_1d3z = tcase_create("NMR PDB-file 1D3Z (several models, hydrogens)");
    tcase_add_test(tc_1d3z, test_1d3z);

    suite_add_tcase(s, tc_basic);
    suite_add_tcase(s, tc_lr_basic);
    suite_add_tcase(s, tc_lr_static);
    suite_add_tcase(s, tc_sr_basic);
    suite_add_tcase(s, tc_lr);
    suite_add_tcase(s, tc_sr);
    suite_add_tcase(s, tc_trimmed);
    suite_add_tcase(s, tc_1d3z);

#if USE_THREADS
    printf("Using pthread\n");
    TCase *tc_pthr = tcase_create("Pthread");
    tcase_add_test(tc_pthr, test_multi_calc);
    suite_add_tcase(s, tc_pthr);
#endif
    return s;
}
