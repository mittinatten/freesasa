#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <check.h>
#include <freesasa.h>
#include <json-c/json_object.h>
#include <json-c/json_object_iterator.h>

#include "tools.h"

extern json_object *
freesasa_node2json(freesasa_node *node, int exclude_type, int options);

static int
compare_nodearea(json_object *obj, const freesasa_nodearea *ref, int is_abs)
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
        ck_assert(float_eq(total, polar + apolar, 1e-10));
        ck_assert(float_eq(total, sc + bb, 1e-10));
        ck_assert(float_eq(total, ref->total, 1e-10));
        ck_assert(float_eq(polar, ref->polar, 1e-10));
        ck_assert(float_eq(apolar, ref->apolar, 1e-10));
        ck_assert(float_eq(sc, ref->side_chain, 1e-10));
        ck_assert(float_eq(bb, ref->main_chain, 1e-10));
    }
    return 1;
}

int test_atom(freesasa_node *node)
{
    ck_assert_ptr_ne(node, NULL);
    json_object *atom = freesasa_node2json(node, FREESASA_NODE_NONE, 0);
    ck_assert_ptr_ne(atom, NULL);

    struct json_object_iterator it = json_object_iter_begin(atom),
                                it_end = json_object_iter_end(atom);
    while (!json_object_iter_equal(&it, &it_end)) {
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

    return 1;
}

int test_residue(freesasa_node *node)
{
    ck_assert_ptr_ne(node, NULL);
    json_object *residue = freesasa_node2json(node, FREESASA_NODE_NONE, 0);
    ck_assert_ptr_ne(residue, NULL);
    const freesasa_nodearea *resarea = freesasa_node_area(node);
    struct json_object_iterator it = json_object_iter_begin(residue),
                                it_end = json_object_iter_end(residue);

    while (!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "name")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "MET");
        } else if (!strcmp(key, "number")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "1");
        } else if (!strcmp(key, "n-atoms")) {
            ck_assert(json_object_is_type(val, json_type_int));
            ck_assert_int_eq(json_object_get_int(val), 8);
        } else if (!strcmp(key, "atoms")) {
            ck_assert(json_object_is_type(val, json_type_array));
            //this is checked further by test_atom
        } else if (!strcmp(key, "area")) {
            ck_assert(compare_nodearea(val, resarea, 1));
        } else if (!strcmp(key, "relative-area")) {
            ck_assert(compare_nodearea(val, resarea, 0));
        } else {
            ck_assert_str_eq(key, "unknown-key");
        }

        json_object_iter_next(&it);
    }

    json_object_put(residue);
    return 1;
}

int test_chain(freesasa_node *node, const freesasa_result *result)
{
    ck_assert_ptr_ne(node, NULL);
    json_object *chain = freesasa_node2json(node, FREESASA_NODE_NONE, 0);
    const freesasa_nodearea *chain_area = freesasa_node_area(node);
    ck_assert_ptr_ne(chain, NULL);
    ck_assert(float_eq(chain_area->total, result->total, 1e-10));

    struct json_object_iterator it = json_object_iter_begin(chain),
                                it_end = json_object_iter_end(chain);
    while (!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "label")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "A");
        } else if (!strcmp(key, "n-residues")) {
            ck_assert(json_object_is_type(val, json_type_int));
            ck_assert_int_eq(json_object_get_int(val), 76);
        } else if (!strcmp(key, "area")) {
            ck_assert(compare_nodearea(val, chain_area, 1));
        } else if (!strcmp(key, "residues")) {
            ck_assert(json_object_is_type(val, json_type_array));
            // the rest is checked in test_residue
        } else {
            ck_assert_str_eq(key, "unknown-key");
        }
        json_object_iter_next(&it);
    }

    json_object_put(chain);
    return 1;
}

int test_structure(freesasa_node *node)
{
    ck_assert_ptr_ne(node, NULL);
    freesasa_nodearea structure_area = {
        .name = "1ubq",
        .total = 4804.0556411417447,
        .polar = 2504.2173023011442,
        .apolar = 2299.838338840601,
        .side_chain = 3689.8982162353718,
        .main_chain = 1114.157424906374};
    json_object *jstruct = freesasa_node2json(node, FREESASA_NODE_NONE, 0);
    ck_assert_ptr_ne(jstruct, NULL);

    struct json_object_iterator it = json_object_iter_begin(jstruct),
                                it_end = json_object_iter_end(jstruct);
    while (!json_object_iter_equal(&it, &it_end)) {
        const char *key = json_object_iter_peek_name(&it);
        json_object *val = json_object_iter_peek_value(&it);
        if (!strcmp(key, "input")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "test");
        } else if (!strcmp(key, "chain-labels")) {
            ck_assert(json_object_is_type(val, json_type_string));
            ck_assert_str_eq(json_object_get_string(val), "A");
        } else if (!strcmp(key, "area")) {
            compare_nodearea(val, &structure_area, 1);
        } else if (!strcmp(key, "model")) {
            ck_assert(json_object_is_type(val, json_type_int));
            // ck_assert_str_eq(json_object_get_string(val), "1");
            // these components are tested in test_chains
        } else if (!strcmp(key, "chains")) {
            ck_assert(json_object_is_type(val, json_type_array));
            // these components are tested in test_chains
        } else {
            ck_assert_str_eq(key, "unknown-key");
        }

        json_object_iter_next(&it);
    }

    json_object_put(jstruct);
    return 1;
}

START_TEST(test_json)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *ubq =
        freesasa_structure_from_pdb(pdb, &freesasa_default_classifier, 0);
    fclose(pdb);
    freesasa_result *result = freesasa_calc_structure(ubq, NULL);
    freesasa_node *tree = freesasa_tree_new();
    freesasa_tree_add_result(tree, result, ubq, "test");
    freesasa_node *result_node = freesasa_node_children(tree);
    freesasa_node *structures = freesasa_node_children(result_node);
    freesasa_node *chains = freesasa_node_children(structures);
    freesasa_node *residues = freesasa_node_children(chains);
    freesasa_node *atoms = freesasa_node_children(residues);

    ck_assert(test_atom(atoms));
    ck_assert(test_residue(residues));
    ck_assert(test_chain(chains, result));
    ck_assert(test_structure(structures));

    freesasa_structure_free(ubq);
    freesasa_result_free(result);
    freesasa_node_free(tree);
}
END_TEST

Suite *json_suite()
{
    Suite *s = suite_create("JSON");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_json);

    suite_add_tcase(s, tc_core);

    return s;
}
