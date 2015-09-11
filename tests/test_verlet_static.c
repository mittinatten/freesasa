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
#include <verlet.c>


#define n_atoms 6
#define n_coord 18
static const double v[n_coord] = {0,0,0, 1,1,1, -1,1,-1, 2,0,-2, 2,2,0, -5,5,5};
static const double r[n_atoms]  = {4,2,2,2,2,2};


START_TEST (test_cell) {
    freesasa_coord *coord = freesasa_coord_new();
    freesasa_coord_append(coord,v,n_atoms);
    double r_max = max_array(r,n_atoms);
    ck_assert(fabs(r_max-4) < 1e-10);
    cell_list *cell_list = cell_list_new(r_max,coord);
    ck_assert(cell_list != NULL);
    ck_assert(fabs(cell_list->d - r_max) < 1e-10);
    ck_assert(fabs(cell_list->x_min - (-1 - r_max/2.)) < 1e-10);
    
    cell_list_free(cell_list);
    freesasa_coord_free(coord);
}
END_TEST

Suite* verlet_static_suite() {
    Suite *s = suite_create("Verlet static functions");

    TCase *tc = tcase_create("Basic");
    tcase_add_test(tc,test_cell);
    
    suite_add_tcase(s, tc);
    
    return s;
}
