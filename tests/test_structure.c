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

sasalib_structure_t *s;

extern int sasalib_set_verbosity(int v);

void setup(void)
{
    s = sasalib_structure_init();
    for (int i = 0; i < N; ++i) {
	sasalib_structure_add_atom(s,an[i],rna[i],rnu[i],cl[i],
				   i,i,i);
    }
}
void teardown(void)
{
    sasalib_structure_free(s);
    s = NULL;
}

START_TEST (test_structure_api)
{

    for (int i = 0; i < N; ++i) {
	ck_assert_str_eq(sasalib_structure_atom_name(s,i),an[i]);
	ck_assert_str_eq(sasalib_structure_atom_res_name(s,i),rna[i]);
	ck_assert_str_eq(sasalib_structure_atom_res_number(s,i),rnu[i]);
	ck_assert_int_eq(sasalib_structure_atom_chain(s,i),cl[i]);
    }
    const sasalib_coord_t *c = sasalib_structure_xyz(s);
    for (int i = 0; i < N; ++i) {
	const double *xyz = sasalib_coord_i(c, i);
	ck_assert(fabs(xyz[0]+xyz[1]+xyz[2]-3*i) < 1e-10);
    }
    ck_assert(sasalib_structure_n(s) == N);
}
END_TEST

START_TEST (test_write_no_pdb)
{
    FILE *null = fopen("/dev/null","w"); // won't work on all platforms

    sasalib_set_verbosity(1);
    ck_assert(sasalib_structure_write_pdb_bfactors(NULL,s,bfactors) == SASALIB_FAIL);
    ck_assert(sasalib_structure_write_pdb_bfactors(null,s,NULL) == SASALIB_FAIL);
    ck_assert(sasalib_structure_write_pdb_bfactors(NULL,s,NULL) == SASALIB_FAIL);
    ck_assert(sasalib_structure_write_pdb_bfactors(null,s,bfactors) == SASALIB_FAIL);
    sasalib_set_verbosity(0);

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
    sasalib_structure_r_def(r,s);
    ck_assert(fabs(r[0]-1.55) < 1e-10);
    ck_assert(fabs(r[1]-2.00) < 1e-10);
    ck_assert(fabs(r[2]-1.40) < 1e-10);
    ck_assert(fabs(r[3]-2.00) < 1e-10);
    ck_assert(fabs(r[4]-2.00) < 1e-10);

    sasalib_structure_r(r,s,a2r);
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
    if (s) sasalib_structure_free(s);
    s = sasalib_structure_init_from_pdb(pdb);
    fclose(pdb);
}

void teardown_1ubq(void)
{
    sasalib_structure_free(s);
    s = NULL;
}

START_TEST (test_structure_1ubq) 
{
    ck_assert(sasalib_structure_n(s) == 602);

    // check at random atom to see that parsing was correct
    ck_assert_str_eq(sasalib_structure_atom_res_name(s,8),"GLN");
    ck_assert_str_eq(sasalib_structure_atom_name(s,8)," N  ");
    ck_assert_str_eq(sasalib_structure_atom_res_number(s,8),"   2");
    ck_assert_int_eq(sasalib_structure_atom_chain(s,8),'A');

    // check coordinates of that random atom
    const sasalib_coord_t *c = sasalib_structure_xyz(s);
    ck_assert(c != NULL);
    const double *x = sasalib_coord_i(c,8);
    ck_assert(x != NULL);
    ck_assert(fabs(x[0]-26.335+x[1]-27.770+x[2]-3.258) < 1e-10);
}
END_TEST

START_TEST (test_write_1ubq) {
    FILE *tf = fopen("tmp/dummy_bfactors.pdb","w+"), 
	*ref = fopen("data/reference_bfactors.pdb","r");
    ck_assert(tf != NULL); 
    ck_assert(ref != NULL);
    const size_t n = sasalib_structure_n(s);
    double *b = (double*)malloc(sizeof(double)*n);
    for (int i = 0; i < n; ++i) b[i] = 1.23;

    ck_assert(sasalib_structure_write_pdb_bfactors(tf,s,b) == SASALIB_SUCCESS);
    
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
    sasalib_set_verbosity(1);
    for (int i = 0; i < 3; ++i) {
	FILE *pdb = fopen(file_names[i],"r");
	ck_assert(pdb != NULL);
	sasalib_structure_t *s = sasalib_structure_init_from_pdb(pdb);
	if (result_null[i]) ck_assert(s == NULL);
	else ck_assert(s != NULL);
	fclose(pdb);
	sasalib_structure_free(s);
    }
    sasalib_set_verbosity(0);
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

