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
#include <check.h>
#include "whole_lib_one_file.c"

int is_identical(const double *l1, const double *l2, int n) {
    for (int i = 0; i < n; ++i) {
        if (l1[i] != l2[i]) return 0;
    }
    return 1;
}

int is_sorted(const double *list,int n)
{
    for (int i = 0; i < n - 1; ++i) if (list[2*i] > list[2*i+1]) return 0;
    return 1;
}

START_TEST (test_sort_arcs) {
    double a_ref[] = {0,1,2,3}, b_ref[] = {-2,0,-1,0,-1,1};
    double a1[4] = {0,1,2,3}, a2[4] = {2,3,0,1};
    double b1[6] = {-2,0,-1,0,-1,1}, b2[6] = {-1,1,-2,0,-1,1};
    sort_arcs(a1,2);
    sort_arcs(a2,2);
    sort_arcs(b1,3);
    sort_arcs(b2,3);
    ck_assert(is_sorted(a1,2));
    ck_assert(is_sorted(a2,2));
    ck_assert(is_sorted(b1,3));
    ck_assert(is_sorted(b2,3));
    ck_assert(is_identical(a_ref,a1,4));
    ck_assert(is_identical(a_ref,a2,4));
    ck_assert(is_identical(b_ref,b1,6));
}
END_TEST

START_TEST (test_exposed_arc_length) 
{
    double a1[4] = {0,0.1*TWOPI,0.9*TWOPI,TWOPI}, a2[4] = {0.9*TWOPI,TWOPI,0,0.1*TWOPI};
    double a3[4] = {0,TWOPI,1,2}, a4[4] = {1,2,0,TWOPI};
    double a5[4] = {0.1*TWOPI,0.2*TWOPI,0.5*TWOPI,0.6*TWOPI};
    double a6[4] = {0.1*TWOPI,0.2*TWOPI,0.5*TWOPI,0.6*TWOPI};
    double a7[4] = {0.1*TWOPI,0.3*TWOPI,0.15*TWOPI,0.2*TWOPI};
    double a8[4] = {0.15*TWOPI,0.2*TWOPI,0.1*TWOPI,0.3*TWOPI};
    double a9[10] = {0.05,0.1, 0.5,0.6, 0,0.15, 0.7,0.8, 0.75,TWOPI};
    ck_assert(fabs(exposed_arc_length(a1,2) - 0.8*TWOPI) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a2,2) - 0.8*TWOPI) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a3,2)) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a4,2)) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a5,2) - 0.8*TWOPI) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a6,2) - 0.8*TWOPI) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a7,2) - 0.8*TWOPI) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a8,2) - 0.8*TWOPI) < 1e-10);
    ck_assert(fabs(exposed_arc_length(a9,5) - 0.45) < 1e-10);
    // can't think of anything more qualitatively different here
}
END_TEST

int main(int argc, char **argv) {
    Suite *s = suite_create("sasa_lr.c static functions");

    TCase *tc = tcase_create("Basic");
    tcase_add_test(tc,test_sort_arcs);
    tcase_add_test(tc,test_exposed_arc_length);
    
    suite_add_tcase(s, tc);
    
    SRunner *sr = srunner_create(s);
    
    srunner_run_all(sr,CK_VERBOSE);

    return (srunner_ntests_failed(sr) == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
