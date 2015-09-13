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
    // FILE without models
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    struct file_interval* it;
    int n = get_models(pdb,&it);
    ck_assert_int_eq(n,0);
    ck_assert(it == NULL);
    fclose(pdb);

    // this file has models
    pdb = fopen(DATADIR "2jo4.pdb","r");
    n = get_models(pdb,&it);
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
}
END_TEST

START_TEST (test_get_chains)
{
    // Test a non PDB file
    FILE *pdb = fopen(DATADIR "err.config", "r");
    struct file_interval whole_file = get_whole_file(pdb), *it = NULL;
    int nc = get_chains(pdb, whole_file, &it, 0);
    fclose(pdb);
    ck_assert_int_eq(nc,0);
    ck_assert_ptr_eq(it,NULL);

    // This file only has one chain
    pdb = fopen(DATADIR "1ubq.pdb", "r");
    whole_file = get_whole_file(pdb);
    nc = get_chains(pdb,whole_file,&it,0);
    fclose(pdb);
    ck_assert_int_eq(nc,1);
    ck_assert_ptr_ne(it,NULL);
    free(it);

    // This file has 4 chains
    pdb = fopen(DATADIR "2jo4.pdb","r");
    int nm = get_models(pdb,&it);
    ck_assert_int_eq(nm,10);
    ck_assert_ptr_ne(it,NULL);
    for (int i = 0; i < nm; ++i) {
        struct file_interval *jt = NULL;
        nc = get_chains(pdb,it[i],&jt,0);
        ck_assert_int_eq(nc,4);
        ck_assert_ptr_ne(jt,NULL);
        for (int j = 1; j < nc; ++j) {
            ck_assert_int_gt(jt[j].begin, jt[j-1].end);
        }
        free(jt);
    }
    free(it);
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
