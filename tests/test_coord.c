#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <check.h>
#include <coord.h>

sasalib_coord_t *coord;

static void setup(void)
{
    coord = sasalib_coord_new();
}
static void teardown(void)
{
    if (coord) free(coord);
}

START_TEST (test_coord)
{
    ck_assert(coord != NULL);
    double xyz[9] = {0,0,0, 1,0,0, 0,1,0};
    sasalib_coord_append(coord,(double*)xyz,3);
    ck_assert(fabs(sasalib_coord_dist(coord,0,2)-1) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,1,2)-sqrt(2)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,1,0)-1) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,0,0)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,0,1)-1) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,0,2)-1) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,1,2)-2) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,1,0)-1) < 1e-10);
    double xyz2[3] = {1,1,1};
    sasalib_coord_set_i(coord,2,(double*)xyz2);
    ck_assert(fabs(sasalib_coord_dist(coord,0,0)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,0,2)-sqrt(3)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,1,2)-sqrt(2)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,0,2)-3) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,1,2)-2) < 1e-10);
    sasalib_coord_set_i_xyz(coord,1,-1,-1,-1);
    ck_assert(fabs(sasalib_coord_dist(coord,1,1)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,1,2)-sqrt(12)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,1,0)-sqrt(3)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,1,2)-12) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist2(coord,0,2)-3) < 1e-10);
    
    /* Corrupts memory...  
    sasalib_coord_set_all(coord,(double*)xyz,3);
    ck_assert(fabs(sasalib_coord_dist(coord,0,0)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,0,1)-1) < 1e-10);
    sasalib_coord_set_all(coord,&xyz[1],2);
    ck_assert(fabs(sasalib_coord_dist(coord,1,2)) < 1e-10);
    ck_assert(fabs(sasalib_coord_dist(coord,0,1)-1) < 1e-10);
    */
    
    double x[2] = {2,2}, y[2] = {1,2}, z[2] = {0,1};
    sasalib_coord_append_xyz(coord,(double*)x,(double*)y,(double*)z,2);
}

END_TEST

Suite* coord_suite()
{
    Suite *s = suite_create("Coord");
    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core,setup,teardown);
    tcase_add_test(tc_core,test_coord);
    suite_add_tcase(s,tc_core);
    return s;
}
