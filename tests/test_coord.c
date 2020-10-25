#include <check.h>
#include <coord.h>
#include <freesasa.h>
#include <freesasa_internal.h>
#include <math.h>
#include <stdlib.h>

#include "tools.h"

coord_t *coord;

static void setup(void)
{
    coord = freesasa_coord_new();
}
static void teardown(void)
{
    if (coord) free(coord);
}

START_TEST(test_coord)
{
    ck_assert_ptr_ne(coord, NULL);
    double xyz[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    ck_assert_int_eq(freesasa_coord_append(coord, (double *)xyz, 3), FREESASA_SUCCESS);
    ck_assert(float_eq(freesasa_coord_dist(coord, 0, 2), 1, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist(coord, 1, 2), sqrt(2), 1e-10));
    ck_assert(float_eq(freesasa_coord_dist(coord, 1, 0), 1, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 0), 0, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 1), 1, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 2), 1, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 1, 2), 2, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 1, 0), 1, 1e-10));
    double xyz2[3] = {1, 1, 1};
    freesasa_coord_set_i(coord, 2, (double *)xyz2);
    ck_assert(float_eq(freesasa_coord_dist(coord, 0, 0), 0, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist(coord, 0, 2), sqrt(3), 1e-10));
    ck_assert(float_eq(freesasa_coord_dist(coord, 1, 2), sqrt(2), 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 2), 3, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 1, 2), 2, 1e-10));
    freesasa_coord_set_i_xyz(coord, 1, -1, -1, -1);
    ck_assert(float_eq(freesasa_coord_dist(coord, 1, 1), 0, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist(coord, 1, 2), sqrt(12), 1e-10));
    ck_assert(float_eq(freesasa_coord_dist(coord, 1, 0), sqrt(3), 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 1, 2), 12, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 2), 3, 1e-10));

    freesasa_coord_set_all(coord, &xyz[3], 2);
    ck_assert(float_eq(freesasa_coord_dist(coord, 0, 0), 0, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 1), 2, 1e-10));

    double x[2] = {2, 2}, y[2] = {1, 2}, z[2] = {0, 1};
    ck_assert_int_eq(freesasa_coord_append_xyz(coord, (double *)x, (double *)y, (double *)z, 2), FREESASA_SUCCESS);
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 2), 2, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 1, 2), 4, 1e-10));
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 3), 6, 1e-10));

    freesasa_coord_set_all_xyz(coord, x, y, z, 2);
    ck_assert(freesasa_coord_n(coord) == 2);
    ck_assert(float_eq(freesasa_coord_dist2(coord, 0, 1), 2, 1e-10));

    freesasa_coord_set_length_i(coord, 1, 6);
    const double *ci = freesasa_coord_i(coord, 1);
    ck_assert(ci != NULL);
    ck_assert(float_eq(ci[0], 4, 1e-10) && float_eq(ci[1], 4, 1e-10) && float_eq(ci[2], 2, 1e-10));

    freesasa_coord_set_length_all(coord, 3);
    ck_assert(float_eq(ci[0], 2, 1e-10) && float_eq(ci[1], 2, 1e-10) && float_eq(ci[2], 1, 1e-10));

    coord_t *c2 = freesasa_coord_new();
    ck_assert(c2 != NULL);
    ck_assert_int_eq(freesasa_coord_append(c2, xyz2, 1), FREESASA_SUCCESS);
    ck_assert(float_eq(freesasa_coord_dist2_12(coord, c2, 1, 0), 2, 1e-10));

    coord_t *c3 = freesasa_coord_clone(c2);
    ck_assert(c3 != NULL);
    ck_assert(float_eq(freesasa_coord_dist2_12(c3, c2, 0, 0), 0, 1e-10));
    freesasa_coord_set_length_all(c2, 10);
    ck_assert_int_eq(freesasa_coord_copy(c3, c2), FREESASA_SUCCESS);
    ck_assert(float_eq(freesasa_coord_dist2_12(c3, c2, 0, 0), 0, 1e-10));
    freesasa_coord_free(c2);
}
END_TEST

START_TEST(test_memerr)
{
    set_fail_after(0);
    static double v[18] = {0, 0, 0, 1, 1, 1, -1, 1, -1, 2, 0, -2, 2, 2, 0, -5, 5, 5};
    struct coord_t coord = {.xyz = v, .n = 6, .is_linked = 0};
    coord_t *coord_dyn = freesasa_coord_new();
    set_fail_after(1);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    void *ptr[] = {freesasa_coord_new(),
                   freesasa_coord_clone(&coord),
                   freesasa_coord_new_linked(v, 1)};
    int ret[] = {freesasa_coord_append(coord_dyn, v, 1),
                 freesasa_coord_append_xyz(coord_dyn, v, v + 1, v + 2, 1)};
    set_fail_after(0);
    for (int i = 0; i < sizeof(ptr) / sizeof(void *); ++i)
        ck_assert_ptr_eq(ptr[i], NULL);
    for (int i = 0; i < sizeof(ret) / sizeof(int); ++i)
        ck_assert_int_eq(ret[i], FREESASA_FAIL);

    freesasa_set_verbosity(FREESASA_V_NORMAL);
    freesasa_coord_free(coord_dyn);
    set_fail_after(0);
}
END_TEST

Suite *coord_suite()
{
    Suite *s = suite_create("Coord");
    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, setup, teardown);
    tcase_add_test(tc_core, test_coord);
    tcase_add_test(tc_core, test_memerr);
    suite_add_tcase(s, tc_core);
    return s;
}
