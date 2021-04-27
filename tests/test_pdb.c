#include <check.h>
#include <math.h>
#include <pdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools.h"

START_TEST(test_pdb_empty_lines)
{
    char buf[80];
    double x[3];

    // check string parsing
    ck_assert_int_eq(freesasa_pdb_get_atom_name(buf, ""), FREESASA_FAIL);
    ck_assert_str_eq(buf, "");
    ck_assert_int_eq(freesasa_pdb_get_res_name(buf, ""), FREESASA_FAIL);
    ck_assert_str_eq(buf, "");
    ck_assert_int_eq(freesasa_pdb_get_res_number(buf, ""), FREESASA_FAIL);
    ck_assert_str_eq(buf, "");
    ck_assert_int_eq(freesasa_pdb_get_symbol(buf, ""), FREESASA_FAIL);
    ck_assert_str_eq(buf, "");

    // check coordinate parsing
    ck_assert_int_eq(freesasa_pdb_get_coord(x, ""), FREESASA_FAIL);

    // check label parsing
    ck_assert_int_eq(freesasa_pdb_get_chain_label(""), '\0');
    ck_assert_int_eq(freesasa_pdb_get_alt_coord_label(""), '\0');

    // check element parsing
    ck_assert_int_eq(freesasa_pdb_ishydrogen(""), FREESASA_FAIL);
}
END_TEST

START_TEST(test_pdb_lines)
{
    char buf[80];
    double x[3], v;
    const char *lines[] = {
        "ATOM    585  C   ARG A  74      41.765  34.829  30.944  0.45 36.22           C",
        "ATOM    573  NH1AARG A  72      34.110  28.437  27.768  1.00 35.02           N  ",
        "HETATM  610  O   HOH A  83      27.707  15.908   4.653  1.00 20.30           O  ",
        "ATOM    573  H   ARG A  72      34.110  28.437  27.768  1.00 35.02           H  ",
        "ATOM    585  C   ARG A  74      41.765  34.829          0.45 36.22           C",
        "HETATM 7673 CD    CD A1978      30.426  14.804  -3.685  1.00 39.42          CD",
        "ATOM    573  H   ARG A  72      34.110  28.437  27.768  1.00 35.02            ",
        "ATOM    573  D   ARG A  72      34.110  28.437  27.768  1.00 35.02           D",
        "ATOM    573  D   ARG A  72      34.110  28.437  27.768  1.00 35.02            ",
    };

    //Atom-name
    freesasa_pdb_get_atom_name(buf, lines[0]);
    ck_assert_str_eq(buf, " C  ");
    freesasa_pdb_get_atom_name(buf, lines[1]);
    ck_assert_str_eq(buf, " NH1");
    freesasa_pdb_get_atom_name(buf, lines[2]);
    ck_assert_str_eq(buf, " O  ");

    //Res-name
    freesasa_pdb_get_res_name(buf, lines[0]);
    ck_assert_str_eq(buf, "ARG");
    freesasa_pdb_get_res_name(buf, lines[2]);
    ck_assert_str_eq(buf, "HOH");

    //Res-number
    freesasa_pdb_get_res_number(buf, lines[0]);
    ck_assert_str_eq(buf, "  74 ");
    freesasa_pdb_get_res_number(buf, lines[1]);
    ck_assert_str_eq(buf, "  72 ");
    freesasa_pdb_get_res_number(buf, lines[2]);
    ck_assert_str_eq(buf, "  83 ");

    //coordinates
    freesasa_pdb_get_coord(x, lines[0]);
    ck_assert(float_eq(x[0], 41.765, 1e-6) &&
              float_eq(x[1], 34.829, 1e-6) &&
              float_eq(x[2], 30.944, 1e-6));
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_pdb_get_coord(x, lines[4]), FREESASA_FAIL);
    freesasa_set_verbosity(FREESASA_V_NORMAL);

    //chain label
    ck_assert_int_eq(freesasa_pdb_get_chain_label(lines[0]), 'A');

    // alt coord labels
    ck_assert_int_eq(freesasa_pdb_get_alt_coord_label(lines[0]), ' ');
    ck_assert_int_eq(freesasa_pdb_get_alt_coord_label(lines[1]), 'A');

    // is hydrogen
    ck_assert(!freesasa_pdb_ishydrogen(lines[0]));
    ck_assert(freesasa_pdb_ishydrogen(lines[3]));
    ck_assert(!freesasa_pdb_ishydrogen(lines[5]));
    ck_assert(freesasa_pdb_ishydrogen(lines[6]));
    ck_assert(freesasa_pdb_ishydrogen(lines[7]));
    ck_assert(freesasa_pdb_ishydrogen(lines[8]));

    // symbol
    ck_assert_int_eq(freesasa_pdb_get_symbol(buf, lines[0]), FREESASA_SUCCESS);
    ck_assert_str_eq(buf, " C");
    ck_assert_int_eq(freesasa_pdb_get_symbol(buf, lines[1]), FREESASA_SUCCESS);
    ck_assert_str_eq(buf, " N");

    // B-factor
    ck_assert_int_eq(freesasa_pdb_get_bfactor(&v, lines[0]), FREESASA_SUCCESS);
    ck_assert(float_eq(v, 36.22, 1e-4));

    // Occupancy
    ck_assert_int_eq(freesasa_pdb_get_occupancy(&v, lines[0]), FREESASA_SUCCESS);
    ck_assert(float_eq(v, 0.45, 1e-4));
}
END_TEST

START_TEST(test_get_models)
{
    // FILE without models
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    struct file_range *it;
    int n = freesasa_pdb_get_models(pdb, &it);
    ck_assert_int_eq(n, 0);
    ck_assert(it == NULL);
    fclose(pdb);

    pdb = fopen(DATADIR "model_mismatch.pdb", "r");
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_pdb_get_models(pdb, &it), FREESASA_FAIL);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    fclose(pdb);

    // this file has models
    pdb = fopen(DATADIR "2jo4.pdb", "r");
    n = freesasa_pdb_get_models(pdb, &it);
    ck_assert_int_eq(n, 10);
    for (int i = 0; i < n; ++i) {
        char *line = NULL;
        size_t len;
        ck_assert_int_gt(it[i].end, it[i].begin);
        fseek(pdb, it[i].begin, 0);
        ck_assert(getline(&line, &len, pdb) > 0);
        // each segment should begin with MODEL
        ck_assert(strncmp(line, "MODEL", 5) == 0);
        while (1) {
            ck_assert(getline(&line, &len, pdb) > 0);
            // there should be only one MODEL per model
            ck_assert(strncmp(line, "MODEL", 5) != 0);
            if (ftell(pdb) >= it[i].end) break;
        }
        // the last line of the segment should be ENDMDL
        ck_assert(strncmp(line, "ENDMDL", 6) == 0);
        free(line);
    }
    free(it);
    fclose(pdb);
}
END_TEST

START_TEST(test_get_chains)
{
    // Test a non PDB file
    FILE *pdb = fopen(DATADIR "err.config", "r");
    struct file_range whole_file = freesasa_whole_file(pdb), *it = NULL;
    int nc = freesasa_pdb_get_chains(pdb, whole_file, &it, 0);
    fclose(pdb);
    ck_assert_int_eq(nc, 0);
    ck_assert_ptr_eq(it, NULL);

    // This file only has one chain
    pdb = fopen(DATADIR "1ubq.pdb", "r");
    whole_file = freesasa_whole_file(pdb);
    nc = freesasa_pdb_get_chains(pdb, whole_file, &it, 0);
    fclose(pdb);
    ck_assert_int_eq(nc, 1);
    ck_assert_ptr_ne(it, NULL);
    free(it);

    // This file has 4 chains
    pdb = fopen(DATADIR "2jo4.pdb", "r");
    int nm = freesasa_pdb_get_models(pdb, &it);
    ck_assert_int_eq(nm, 10);
    ck_assert_ptr_ne(it, NULL);
    for (int i = 0; i < nm; ++i) {
        struct file_range *jt = NULL;
        nc = freesasa_pdb_get_chains(pdb, it[i], &jt, 0);
        ck_assert_int_eq(nc, 4);
        ck_assert_ptr_ne(jt, NULL);
        for (int j = 1; j < nc; ++j) {
            ck_assert_int_ge(jt[j].begin, jt[j - 1].end);
        }
        free(jt);
    }
    free(it);
}
END_TEST

extern TCase *test_pdb_static();

Suite *pdb_suite()
{
    Suite *s = suite_create("PDB-parser");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_pdb_empty_lines);
    tcase_add_test(tc_core, test_pdb_lines);
    tcase_add_test(tc_core, test_get_models);
    tcase_add_test(tc_core, test_get_chains);

    TCase *tc_static = test_pdb_static();

    suite_add_tcase(s, tc_core);
    suite_add_tcase(s, tc_static);

    return s;
}
