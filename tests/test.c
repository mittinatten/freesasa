/*
  Copyright Simon Mitternacht 2013-2014.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
//#include <stdio.h>
#include <check.h>

extern int test_sasa_basic();
extern Suite* pdb_suite();

int main(int argc, char **argv) {
    int n_err = 0;
    
    n_err += test_sasa_basic();

    Suite *s = pdb_suite();
    SRunner *sr = srunner_create(s);
    srunner_run_all(sr,CK_NORMAL);
    n_err += srunner_ntests_failed(sr);
    //printf("Total: there were %d errors.\n",n_err);

    return (n_err == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
