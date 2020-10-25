#include <check.h>
#include <freesasa_internal.h>
#include <nb.h>

#include "tools.h"

const double v[18] = {0, 0, 0, 1, 1, 1, -1, 1, -1, 2, 0, -2, 2, 2, 0, -5, 5, 5};
const double r[6] = {4, 2, 2, 2, 2, 2};

START_TEST(test_nb)
{
    coord_t *coord = freesasa_coord_new();
    nb_list *nb;
    freesasa_coord_append(coord, v, 6);
    ck_assert_ptr_eq(freesasa_nb_new(NULL, NULL), NULL);
    ck_assert_ptr_eq(freesasa_nb_new(NULL, r), NULL);
    ck_assert_ptr_eq(freesasa_nb_new(coord, NULL), NULL);

    nb = freesasa_nb_new(coord, r);
    ck_assert(nb != NULL);
    ck_assert(freesasa_nb_contact(nb, 0, 1));
    ck_assert(freesasa_nb_contact(nb, 1, 0));
    ck_assert(freesasa_nb_contact(nb, 0, 5) == 0);
    freesasa_nb_free(nb);
    freesasa_coord_free(coord);
}
END_TEST

START_TEST(test_memerr)
{
    freesasa_set_verbosity(FREESASA_V_SILENT);
    double v[18] = {0, 0, 0, 1, 1, 1, -1, 1, -1, 2, 0, -2, 2, 2, 0, -5, 5, 5};
    struct coord_t coord = {.xyz = v, .n = 6, .is_linked = 0};
    const double r[6] = {4, 2, 2, 2, 2, 2};

    for (int i = 1; i < 30; ++i) {
        set_fail_after(i);
        void *ptr = freesasa_nb_new(&coord, r);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

extern TCase *test_nb_static();

Suite *nb_suite()
{
    Suite *s = suite_create("Neighbor lists");

    TCase *tc_nb = tcase_create("Basic");
    tcase_add_test(tc_nb, test_nb);
    tcase_add_test(tc_nb, test_memerr);

    TCase *tc_static = test_nb_static();

    suite_add_tcase(s, tc_nb);
    suite_add_tcase(s, tc_static);

    return s;
}
