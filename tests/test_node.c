#include <check.h>
#include <freesasa.h>
#include <freesasa_internal.h>

#include "tools.h"

static void
test_tree(freesasa_node *structure,
          const freesasa_result *result)
{
    freesasa_node *next, *chain, *residue, *atom;

    const freesasa_nodearea *area;
    ck_assert_ptr_ne((chain = freesasa_node_children(structure)), NULL);
    ck_assert_ptr_ne((residue = freesasa_node_children(chain)), NULL);
    ck_assert_ptr_ne((atom = freesasa_node_children(residue)), NULL);

    ck_assert_int_eq(freesasa_node_type(structure), FREESASA_NODE_STRUCTURE);
    ck_assert_int_eq(freesasa_node_type(chain), FREESASA_NODE_CHAIN);
    ck_assert_int_eq(freesasa_node_type(residue), FREESASA_NODE_RESIDUE);
    ck_assert_int_eq(freesasa_node_type(atom), FREESASA_NODE_ATOM);

    ck_assert_str_eq(freesasa_node_name(structure), "A");
    ck_assert_str_eq(freesasa_node_name(chain), "A");
    ck_assert_str_eq(freesasa_node_name(residue), "MET");
    ck_assert_str_eq(freesasa_node_name(atom), " N  ");

    ck_assert_int_eq(freesasa_node_structure_n_chains(structure), 1);
    ck_assert_int_eq(freesasa_node_structure_n_atoms(structure), 602);
    ck_assert_str_eq(freesasa_node_structure_chain_labels(structure), "A");

    ck_assert_int_eq(freesasa_node_chain_n_residues(chain), 76);

    // iterate
    next = freesasa_node_next(structure);
    ck_assert_ptr_eq(next, NULL);
    next = freesasa_node_next(chain);
    ck_assert_ptr_eq(next, NULL);
    next = freesasa_node_next(residue);
    ck_assert_ptr_ne(next, NULL);
    ck_assert_str_eq(freesasa_node_name(next), "GLN");
    next = freesasa_node_next(atom);
    ck_assert_str_eq(freesasa_node_name(next), " CA ");

    ck_assert_ptr_ne((area = freesasa_node_area(structure)), NULL);
    ck_assert(float_eq(result->total, area->total, 1e-10));
    ck_assert_ptr_ne((area = freesasa_node_area(chain)), NULL);
    ck_assert(float_eq(result->total, area->total, 1e-10));
    next = freesasa_node_next(residue);
    ck_assert_ptr_ne((area = freesasa_node_area(next)), NULL);
    ck_assert_ptr_ne((area = freesasa_node_area(atom)), NULL);
    ck_assert(float_eq(result->sasa[0], area->total, 1e-10));
    next = freesasa_node_next(atom);
    ck_assert_ptr_ne((area = freesasa_node_area(next)), NULL);
    ck_assert(float_eq(result->sasa[1], area->total, 1e-10));
}

START_TEST(test_result_node)
{
    FILE *file = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *structure = freesasa_structure_from_pdb(file, NULL, 0);
    freesasa_result *result = freesasa_calc_structure(structure, NULL);
    freesasa_node *tree = freesasa_tree_new(), *tree2 = freesasa_tree_new();
    freesasa_node *rn;

    ck_assert_ptr_ne(tree, NULL);
    ck_assert_ptr_ne(tree2, NULL);

    ck_assert_int_eq(freesasa_tree_add_result(tree, result, structure, "test"), FREESASA_SUCCESS);
    ck_assert_int_eq(freesasa_tree_add_result(tree2, result, structure, "test2"), FREESASA_SUCCESS);

    ck_assert_ptr_ne(tree, NULL);
    ck_assert_ptr_eq(freesasa_node_parent(tree), NULL);
    ck_assert_int_eq(freesasa_node_type(tree), FREESASA_NODE_ROOT);
    ck_assert_ptr_eq(freesasa_node_name(tree), NULL);
    ck_assert_ptr_eq(freesasa_node_next(tree), NULL);
    ck_assert_ptr_ne((rn = freesasa_node_children(tree)), NULL);
    ck_assert_str_eq(freesasa_node_name(rn), "test");
    ck_assert_int_eq(freesasa_node_type(rn), FREESASA_NODE_RESULT);
    ck_assert_ptr_ne((rn = freesasa_node_children(rn)), NULL);
    test_tree(rn, result);

    freesasa_tree_join(tree, &tree2);
    ck_assert(tree2 == NULL);
    rn = freesasa_node_children(tree); // result in tree
    ck_assert_ptr_ne(rn, NULL);
    rn = freesasa_node_children(rn); // structure in tree
    test_tree(rn, result);
    rn = freesasa_node_children(tree);
    rn = freesasa_node_next(rn); // result in tree2
    ck_assert_ptr_ne(rn, NULL);
    rn = freesasa_node_children(rn); // structure in tree2
    ck_assert_ptr_ne(rn, NULL);
    test_tree(rn, result);

    freesasa_node_free(tree);
    freesasa_structure_free(structure);
    freesasa_result_free(result);
}
END_TEST

START_TEST(test_memerr)
{
    FILE *file = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *structure = freesasa_structure_from_pdb(file, NULL, 0);
    freesasa_result *result = freesasa_calc_structure(structure, NULL);
    freesasa_node *rn;
    freesasa_set_verbosity(FREESASA_V_SILENT);
    rn = freesasa_tree_new();
    for (int i = 1; i < 200; ++i) {
        int ret;
        set_fail_after(i);
        ret = freesasa_tree_add_result(rn, result, structure, "test");
        set_fail_after(0);
        ck_assert_int_eq(ret, FREESASA_FAIL);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    freesasa_structure_free(structure);
    freesasa_result_free(result);
    fclose(file);
}
END_TEST

Suite *result_node_suite()
{
    Suite *s = suite_create("Result-node");

    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_result_node);
    tcase_add_test(tc_core, test_memerr);

    suite_add_tcase(s, tc_core);

    return s;
}
