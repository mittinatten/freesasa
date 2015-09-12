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

START_TEST (test_get_models) {
    // this file has models
    FILE *pdb = fopen(DATADIR "2jo4.pdb","r");
    struct file_interval* it;
    int n = get_models(pdb,&it);
    ck_assert_int_eq(n,10);
    for (int i = 0; i < n; ++i) {
        char *line = NULL;
        size_t len;
        ck_assert_int_gt(it[i].end,it[i].begin);
        fseek(pdb,it[i].begin,0);
        getline(&line,&len,pdb);
        // each segment should begin with MODEL
        ck_assert(strncmp(line,"MODEL",5) == 0);
        while(1) {
            getline(&line,&len,pdb);
            // there should be only one MODEL per model
            ck_assert(strncmp(line,"MODEL",5) != 0);
            if (ftell(pdb) >= it[i].end) break;
        }
        // the last line of the segment should be ENDMDL
        ck_assert(strncmp(line,"ENDMDL",6) == 0);
        free(line);
    }
    free(it);
    fclose(pdb);
    pdb = fopen(DATADIR "1ubq.pdb","r");
    n = get_models(pdb,&it);
    ck_assert_int_eq(n,0);
    ck_assert(it == NULL);
    fclose(pdb);
}
END_TEST

int main(int argc, char **argv) {
    Suite *s = suite_create("structure.c static functions");

    TCase *tc = tcase_create("Basic");
    tcase_add_test(tc,test_get_models);
    
    suite_add_tcase(s, tc);
    
    SRunner *sr = srunner_create(s);
    
    srunner_run_all(sr,CK_VERBOSE);

    return (srunner_ntests_failed(sr) == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
