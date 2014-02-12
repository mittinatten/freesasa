#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <check.h>
#include <structure.h>

START_TEST (test_structure_api)
{
    sasalib_structure_t *s = sasalib_structure_init();
#define NRES 5
    const char an[NRES][5] =  {" C  "," CA "," O  "," CB "," SD "};
    const char rna[NRES][4] = {"MET", "MET", "MET", "MET", "MET"};
    const char rnu[NRES][5] = {"   1","   1","   1","   1","   1"};
    const char cl[NRES]   = {'A','A','A','A','A'};
    for (int i = 0; i < NRES; ++i) {
	sasalib_structure_add_atom(s,an[i],rna[i],rnu[i],cl[i],
				   i,i,i);
    }
    for (int i = 0; i < NRES; ++i) {
	ck_assert_str_eq(sasalib_structure_atom_name(s,i),an[i]);
	ck_assert_str_eq(sasalib_structure_atom_res_name(s,i),rna[i]);
	ck_assert_str_eq(sasalib_structure_atom_res_number(s,i),rnu[i]);
	ck_assert_int_eq(sasalib_structure_atom_chain(s,i),cl[i]);
    }
    const sasalib_coord_t *c = sasalib_structure_xyz(s);
    for (int i = 0; i < NRES; ++i) {
	const double *xyz = sasalib_coord_i(c, i);
	ck_assert(fabs(xyz[0]+xyz[1]+xyz[2]-3*i) < 1e-10);
    }
    
    sasalib_structure_free(s);
}
END_TEST

Suite* structure_suite() {
    Suite *s = suite_create("Structure");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_structure_api);
    suite_add_tcase(s, tc_core);

    return s;
}
