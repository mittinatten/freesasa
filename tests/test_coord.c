/*
  Copyright Simon Mitternacht 2013-2016.

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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <check.h>
#include <freesasa.h>
#include <coord.h>

coord_t *coord;

static void setup(void)
{
    coord = freesasa_coord_new();
}
static void teardown(void)
{
    if (coord) free(coord);
}

START_TEST (test_coord)
{
    ck_assert_ptr_ne(coord,NULL);
    double xyz[9] = {0,0,0, 1,0,0, 0,1,0};
    ck_assert_int_eq(freesasa_coord_append(coord,(double*)xyz,3),FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_coord_dist(coord,0,2)-1) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist(coord,1,2)-sqrt(2)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist(coord,1,0)-1) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,0)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,1)-1) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,2)-1) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,1,2)-2) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,1,0)-1) < 1e-10);
    double xyz2[3] = {1,1,1};
    freesasa_coord_set_i(coord,2,(double*)xyz2);
    ck_assert(fabs(freesasa_coord_dist(coord,0,0)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist(coord,0,2)-sqrt(3)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist(coord,1,2)-sqrt(2)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,2)-3) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,1,2)-2) < 1e-10);
    freesasa_coord_set_i_xyz(coord,1,-1,-1,-1);
    ck_assert(fabs(freesasa_coord_dist(coord,1,1)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist(coord,1,2)-sqrt(12)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist(coord,1,0)-sqrt(3)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,1,2)-12) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,2)-3) < 1e-10);

    freesasa_coord_set_all(coord,&xyz[3],2);
    ck_assert(fabs(freesasa_coord_dist(coord,0,0)) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,1)-2) < 1e-10);

    double x[2] = {2,2}, y[2] = {1,2}, z[2] = {0,1};
    ck_assert_int_eq(freesasa_coord_append_xyz(coord,(double*)x,(double*)y,(double*)z,2),FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,2)-2) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,1,2)-4) < 1e-10);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,3)-6) < 1e-10);

    freesasa_coord_set_all_xyz(coord,x,y,z,2);
    ck_assert(freesasa_coord_n(coord) == 2);
    ck_assert(fabs(freesasa_coord_dist2(coord,0,1)-2) < 1e-10);

    freesasa_coord_set_length_i(coord,1,6);
    const double *ci = freesasa_coord_i(coord,1);
    ck_assert(ci != NULL);
    ck_assert(ci[0] == 4 && ci[1] == 4 && ci[2] == 2);

    freesasa_coord_set_length_all(coord,3);
    ck_assert(ci[0] == 2 && ci[1] == 2 && ci[2] == 1);

    coord_t *c2 = freesasa_coord_new();
    ck_assert(c2 != NULL);
    ck_assert_int_eq(freesasa_coord_append(c2,xyz2,1),FREESASA_SUCCESS);
    ck_assert(fabs(freesasa_coord_dist2_12(coord,c2,0,0) - 2));
    freesasa_coord_free(c2);
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
