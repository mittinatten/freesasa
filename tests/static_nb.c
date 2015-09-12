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
#include <coord.c>
#include <nb.c>


const int n_atoms = 6;
const int n_coord = 18;
static const double v[] = {0,0,0, 1,1,1, -1,1,-1, 2,0,-2, 2,2,0, -5,5,5};
static const double r[]  = {4,2,2,2,2,2};


START_TEST (test_cell) {
    double r_max;
    cell_list *c;
    freesasa_coord *coord = freesasa_coord_new();
    freesasa_coord_append(coord,v,n_atoms);
    r_max = max_array(r,n_atoms);
    ck_assert(fabs(r_max-4) < 1e-10);
    c = cell_list_new(r_max,coord);
    ck_assert(c != NULL);
    ck_assert(c->cell != NULL);
    ck_assert(fabs(c->d - r_max) < 1e-10);
    // check bounds
    ck_assert(fabs(c->x_min < -5));
    ck_assert(fabs(c->x_max > 2));
    ck_assert(fabs(c->y_min < 0));
    ck_assert(fabs(c->y_max > 5));
    ck_assert(fabs(c->z_min < -2));
    ck_assert(fabs(c->z_max > 5));
    // check number of cells
    ck_assert(c->nx*c->d >= 7);
    ck_assert(c->nx <= ceil(7/r_max)+1);
    ck_assert(c->ny*c->d >= 5);
    ck_assert(c->ny <= ceil(5/r_max)+1);
    ck_assert(c->nz*c->d >= 7);
    ck_assert(c->nz <= ceil(7/r_max)+1);
    ck_assert_int_eq(c->n, c->nx*c->ny*c->nz);
    // check the individual cells
    int na = 0;
    ck_assert_int_eq(c->cell[0].n_nb,8);
    ck_assert_int_eq(c->cell[c->n-1].n_nb,1);
    for (int i = 0; i < c->n; ++i) {
        cell ci = c->cell[i];
        ck_assert(ci.n_atoms >= 0);
        if (ci.n_atoms > 0) ck_assert(ci.atom != NULL);
        ck_assert_int_ge(ci.n_nb, 1); 
        ck_assert_int_le(ci.n_nb, 17);
        na += ci.n_atoms;
    }
    ck_assert_int_eq(na,n_atoms);
    cell_list_free(c);
    freesasa_coord_free(coord);
}
END_TEST

int main(int argc, char **argv) {
    Suite *s = suite_create("nb.c static functions");

    TCase *tc = tcase_create("Basic");
    tcase_add_test(tc,test_cell);
    
    suite_add_tcase(s, tc);
    
    SRunner *sr = srunner_create(s);
    
    srunner_run_all(sr,CK_VERBOSE);

    return (srunner_ntests_failed(sr) == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
