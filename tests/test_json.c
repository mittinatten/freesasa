#include <check.h>
#include <freesasa_json.h>
#include <json-c/json_object_iterator.h>
#include <json-c/json_util.h>
#include "tools.h"

static int
compare_rs(json_object *obj, freesasa_residue_sasa *ref, int is_abs)
{
    struct json_object_iterator it = json_object_iter_begin(obj),
        it_end = json_object_iter_end(obj);
    double total, polar, apolar, bb, sc;
    while (!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "total")) {
            ck_assert(json_object_is_type(val, json_type_double));
            total = json_object_get_double(val);
        } else if (!strcmp(key, "polar")) {
            ck_assert(json_object_is_type(val, json_type_double));
            polar = json_object_get_double(val);
        } else if (!strcmp(key, "apolar")) {
            ck_assert(json_object_is_type(val, json_type_double));
            apolar = json_object_get_double(val);
        } else if (!strcmp(key, "main-chain")) {
            ck_assert(json_object_is_type(val, json_type_double));
            bb = json_object_get_double(val);
        } else if (!strcmp(key, "side-chain")) {
            ck_assert(json_object_is_type(val, json_type_double));
            sc = json_object_get_double(val);
        } else {
            ck_assert(0);
        }
        json_object_iter_next(&it);
    }
    ck_assert(total > 0);
    if (is_abs) {
        ck_assert(float_eq(total, polar+apolar, 1e-10));
        ck_assert(float_eq(total, sc+bb, 1e-10));
        ck_assert(float_eq(total, ref->total, 1e-10));
        ck_assert(float_eq(polar, ref->polar, 1e-10));
        ck_assert(float_eq(apolar, ref->apolar, 1e-10));
        ck_assert(float_eq(sc, ref->side_chain, 1e-10));
        ck_assert(float_eq(bb, ref->main_chain, 1e-10));
    }
    return 1;
}

START_TEST (test_atom)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *ubq =
        freesasa_structure_from_pdb(pdb, &freesasa_default_classifier, 0);
    fclose(pdb);
    freesasa_result *result = freesasa_calc_structure(ubq, NULL);
    json_object *atom = freesasa_json_atom(result, ubq, &freesasa_default_rsa, 0);

    ck_assert_ptr_ne(atom, NULL);
    
    struct json_object_iterator it = json_object_iter_begin(atom),
        it_end = json_object_iter_end(atom);
    while(!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "name")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "N");
        } else if (!strcmp(key, "area")) {
            ck_assert(json_object_is_type(val, json_type_double));
            ck_assert(json_object_get_double(val) > 0);
        } else if (!strcmp(key, "is-polar")) {
            ck_assert(json_object_is_type(val, json_type_boolean));
            ck_assert(json_object_get_boolean(val));
        } else if (!strcmp(key, "is-main-chain")) {
            ck_assert(json_object_is_type(val, json_type_boolean));
            ck_assert(json_object_get_boolean(val));
        } else if (!strcmp(key, "radius")) {
            ck_assert(json_object_is_type(val, json_type_double));
            ck_assert(json_object_get_double(val) > 0);
        } else {
            ck_assert(0);
        }
        json_object_iter_next(&it);
    }
    json_object_put(atom);
    freesasa_structure_free(ubq);
} 
END_TEST

START_TEST (test_residue)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *ubq =
        freesasa_structure_from_pdb(pdb, &freesasa_default_classifier, 0);
    fclose(pdb);
    freesasa_result *result = freesasa_calc_structure(ubq, NULL);
    freesasa_residue_sasa chrs = {"A", 0, 0, 0, 0, 0};
    json_object *residue = freesasa_json_residue(result, ubq, &freesasa_default_rsa,
                                                 0, &chrs);

    struct json_object_iterator it = json_object_iter_begin(residue),
        it_end = json_object_iter_end(residue);
    while(!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "name")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "MET");
        } else if (!strcmp(key, "number")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "1");
        } else if (!strcmp(key, "n_atoms")) {
            ck_assert(json_object_is_type(val, json_type_int));
            ck_assert_int_eq(json_object_get_int(val), 8);
        } else if (!strcmp(key, "atoms")) {
             ck_assert(json_object_is_type(val, json_type_array));
            //this is checked further by test_atom
        } else if (!strcmp(key, "abs")) {
            ck_assert(compare_rs(val, &chrs, 1));
        } else if (!strcmp(key, "rel")) {
            ck_assert(compare_rs(val, &chrs, 0));
        } else {
            ck_assert(0);
        }

        json_object_iter_next(&it);
    }

    json_object_put(residue);
    freesasa_structure_free(ubq);
}
END_TEST

START_TEST (test_chain)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *ubq =
        freesasa_structure_from_pdb(pdb, &freesasa_default_classifier, 0);
    fclose(pdb);
    freesasa_result *result = freesasa_calc_structure(ubq, NULL);
    freesasa_residue_sasa structure_rs = {"1ubq", 0, 0, 0, 0, 0};

    json_object *chain = freesasa_json_chain(result, ubq, &freesasa_default_rsa,
                                             'A', &structure_rs);
    
    ck_assert(float_eq(structure_rs.total, result->total, 1e-10));

    struct json_object_iterator it = json_object_iter_begin(chain),
        it_end = json_object_iter_end(chain);
    while(!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "label")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "A");
        } else if (!strcmp(key, "n_residues")) {
            ck_assert(json_object_is_type(val, json_type_int));
            ck_assert_int_eq(json_object_get_int(val), 76);
        } else if (!strcmp(key, "abs")) {
            ck_assert(compare_rs(val, &structure_rs, 1));
        } else if (!strcmp(key, "residues")) {
            ck_assert(json_object_is_type(val, json_type_array));
            // the rest is checked in test_residue
        } else {
            ck_assert(0);
        }
        json_object_iter_next(&it);
    }

    json_object_put(chain);
    freesasa_structure_free(ubq);
    freesasa_result_free(result);
}
END_TEST

START_TEST (test_structure)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *ubq =
        freesasa_structure_from_pdb(pdb, &freesasa_default_classifier, 0);
    fclose(pdb);
    freesasa_result *result = freesasa_calc_structure(ubq, NULL);
    freesasa_residue_sasa structure_rs = {
        .name = "1ubq",
        .total = 4804.0556411417447,
        .polar = 2504.2173023011442,
        .apolar = 2299.838338840601,
        .side_chain = 3689.8982162353718,
        .main_chain = 1114.157424906374
    };
        
    json_object *jstruct = freesasa_json_result(result, ubq, &freesasa_default_rsa, "1ubq");

    struct json_object_iterator it = json_object_iter_begin(jstruct),
        it_end = json_object_iter_end(jstruct);
    while(!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "name")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "1ubq");
        } else if (!strcmp(key, "n_chains")) {
            ck_assert(json_object_is_type(val, json_type_int));
            ck_assert_int_eq(json_object_get_int(val), 1);
        } else if (!strcmp(key, "abs")) {
            compare_rs(val, &structure_rs, 1);
        } else if (!strcmp(key, "chains")) {
            ck_assert(json_object_is_type(val, json_type_array));
            // these components are tested in test_chains
        } else {
            ck_assert(0);
        }

        json_object_iter_next(&it);
    }
    
    json_object_put(jstruct);
    freesasa_structure_free(ubq);
    freesasa_result_free(result);
}
END_TEST

Suite* json_suite() {
    Suite *s = suite_create("JSON");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_atom);
    tcase_add_test(tc_core, test_residue);
    tcase_add_test(tc_core, test_chain);
    tcase_add_test(tc_core, test_structure);

    suite_add_tcase(s, tc_core);

    return s;
}
