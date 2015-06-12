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
#include <adjacency.h>
#include <check.h>

const double v[18] = {0,0,0, 1,1,1, -1,1,-1, 2,0,-2, 2,2,0, -5,5,5};
const double r[6]  = {4,2,2,2,2,2};

START_TEST (test_adjacency) {
    freesasa_coord *coord = freesasa_coord_new();
    freesasa_coord_append(coord,v,6);
    freesasa_adjacency *adj = freesasa_adjacency_new(coord,r);
    ck_assert(adj != NULL);
    ck_assert(freesasa_adjacency_contact(adj,0,1));
    ck_assert(freesasa_adjacency_contact(adj,1,0));
    ck_assert(freesasa_adjacency_contact(adj,0,5) == 0);
    freesasa_adjacency_free(adj);
    freesasa_coord_free(coord);
}
END_TEST

Suite* adjacency_suite() {
    Suite *s = suite_create("Adjacency");

    TCase *tc_adjacency = tcase_create("Basic");
    tcase_add_test(tc_adjacency,test_adjacency);
    
    suite_add_tcase(s, tc_adjacency);
    
    return s;
}
