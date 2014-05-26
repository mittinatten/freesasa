/*
  Copyright Simon Mitternacht 2013-2014.

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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <check.h>
#include <structure.h>
#include <pdb.h>

#define N 5
const char an[N][PDB_ATOM_NAME_STRL+1] =  {" C  "," CA "," O  "," CB "," SD "};
const char rna[N][PDB_ATOM_RES_NAME_STRL+1] = {"MET", "MET", "MET", "MET", "MET"};
const char rnu[N][PDB_ATOM_RES_NUMBER_STRL+1] = {"   1","   1","   1","   1","   1"};
const char cl[N] = {'A','A','A','A','A'};
const double bfactors[N] = {1., 1., 1., 1., 1.};

freesasa_structure_t *s;

extern int freesasa_set_verbosity(int v);

static void setup(void)
{
    s = freesasa_structure_init();
    for (int i = 0; i < N; ++i) {
	freesasa_structure_add_atom(s,an[i],rna[i],rnu[i],cl[i],
				   i,i,i);
    }
}
static void teardown(void)
{
    freesasa_structure_free(s);
    s = NULL;
}

START_TEST (test_structure_api)
{

    for (int i = 0; i < N; ++i) {
	ck_assert_str_eq(freesasa_structure_atom_name(s,i),an[i]);
	ck_assert_str_eq(freesasa_structure_atom_res_name(s,i),rna[i]);
	ck_assert_str_eq(freesasa_structure_atom_res_number(s,i),rnu[i]);
	ck_assert_int_eq(freesasa_structure_atom_chain(s,i),cl[i]);
    }
    const freesasa_coord_t *c = freesasa_structure_xyz(s);
    for (int i = 0; i < N; ++i) {
	const double *xyz = freesasa_coord_i(c, i);
	ck_assert(fabs(xyz[0]+xyz[1]+xyz[2]-3*i) < 1e-10);
    }
    ck_assert(freesasa_structure_n(s) == N);
}
END_TEST

START_TEST (test_write_no_pdb)
{
    FILE *null = fopen("/dev/null","w"); // won't work on all platforms

    freesasa_set_verbosity(1);
    ck_assert(freesasa_structure_write_pdb_bfactors(NULL,s,bfactors) == FREESASA_FAIL);
    ck_assert(freesasa_structure_write_pdb_bfactors(null,s,NULL) == FREESASA_FAIL);
    ck_assert(freesasa_structure_write_pdb_bfactors(NULL,s,NULL) == FREESASA_FAIL);
    ck_assert(freesasa_structure_write_pdb_bfactors(null,s,bfactors) == FREESASA_FAIL);
    freesasa_set_verbosity(0);

    fclose(null);
}
END_TEST

double a2r(const char *rn, const char *am)
{
    return 1.0;
}

/* this will only test that reasonable values are obtained, the tests
   for classify.h will do more exhaustive classification analysis */
START_TEST (test_radii) {
    double r[N];
    freesasa_structure_r_def(r,s);
    ck_assert(fabs(r[0]-1.55) < 1e-10);
    ck_assert(fabs(r[1]-2.00) < 1e-10);
    ck_assert(fabs(r[2]-1.40) < 1e-10);
    ck_assert(fabs(r[3]-2.00) < 1e-10);
    ck_assert(fabs(r[4]-2.00) < 1e-10);

    freesasa_structure_r(r,s,a2r);
    for (int i = 0; i < N; ++i) ck_assert(fabs(r[i]-1.0) < 1e-10);
}
END_TEST

void setup_1ubq(void)
{
    errno = 0;
    FILE *pdb = fopen("data/1ubq.pdb","r");
    if (pdb == NULL) {
    	fprintf(stderr,"error reading PDB-file for test. "
		"(Tests must be run from directory test/): %s\n",
		strerror(errno));
    }
    if (s) freesasa_structure_free(s);
    s = freesasa_structure_init_from_pdb(pdb);
    fclose(pdb);
}

void teardown_1ubq(void)
{
    freesasa_structure_free(s);
    s = NULL;
}

START_TEST (test_structure_1ubq) 
{
    ck_assert(freesasa_structure_n(s) == 602);

    // check at random atom to see that parsing was correct
    ck_assert_str_eq(freesasa_structure_atom_res_name(s,8),"GLN");
    ck_assert_str_eq(freesasa_structure_atom_name(s,8)," N  ");
    ck_assert_str_eq(freesasa_structure_atom_res_number(s,8),"   2");
    ck_assert_int_eq(freesasa_structure_atom_chain(s,8),'A');

    // check coordinates of that random atom
    const freesasa_coord_t *c = freesasa_structure_xyz(s);
    ck_assert(c != NULL);
    const double *x = freesasa_coord_i(c,8);
    ck_assert(x != NULL);
    ck_assert(fabs(x[0]-26.335+x[1]-27.770+x[2]-3.258) < 1e-10);
}
END_TEST

START_TEST (test_write_1ubq) {
    FILE *tf = fopen("tmp/dummy_bfactors.pdb","w+"), 
	*ref = fopen("data/reference_bfactors.pdb","r");
    ck_assert(tf != NULL); 
    ck_assert(ref != NULL);
    const size_t n = freesasa_structure_n(s);
    double *b = (double*)malloc(sizeof(double)*n);
    for (int i = 0; i < n; ++i) b[i] = 1.23;

    ck_assert(freesasa_structure_write_pdb_bfactors(tf,s,b) == FREESASA_SUCCESS);
    
    rewind(tf);

    //check that output matches refernce file
    size_t bufsize = 100;
    char *buf_tf = malloc(bufsize), *buf_ref = malloc(bufsize);
    while(getline(&buf_tf,&bufsize,tf) > 0 && getline(&buf_ref,&bufsize,ref) > 0) { 
	ck_assert_str_eq(buf_ref,buf_tf);
    }
    free(buf_tf);
    free(buf_ref);
    fclose(ref);
    fclose(tf);
}
END_TEST

START_TEST (test_pdb)
{
    const char *file_names[] = {"data/alt_model_twochain.pdb",
				"data/empty.pdb",
				"data/empty_model.pdb"};
    const int result_null[] = {0,1,1};
    freesasa_set_verbosity(1);
    for (int i = 0; i < 3; ++i) {
	FILE *pdb = fopen(file_names[i],"r");
	ck_assert(pdb != NULL);
	freesasa_structure_t *s = freesasa_structure_init_from_pdb(pdb);
	if (result_null[i]) ck_assert(s == NULL);
	else ck_assert(s != NULL);
	fclose(pdb);
	freesasa_structure_free(s);
    }
    freesasa_set_verbosity(0);
}
END_TEST

Suite* structure_suite() {
    Suite *s = suite_create("Structure");
    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core,setup,teardown);
    tcase_add_test(tc_core, test_structure_api);
    tcase_add_test(tc_core, test_radii);
    tcase_add_test(tc_core, test_write_no_pdb);

    TCase *tc_pdb = tcase_create("PDB");
    tcase_add_test(tc_pdb,test_pdb);

    TCase *tc_1ubq = tcase_create("1UBQ");
    tcase_add_checked_fixture(tc_1ubq,setup_1ubq,teardown_1ubq);
    tcase_add_test(tc_1ubq,test_structure_1ubq);
    tcase_add_test(tc_1ubq,test_write_1ubq);
    
    suite_add_tcase(s, tc_core);
    suite_add_tcase(s, tc_pdb);
    suite_add_tcase(s, tc_1ubq);

    return s;
}

