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

freesasa_structure *s;

static void setup(void)
{
    s = freesasa_structure_new();
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
    char buf[20];
    sprintf(buf,"%c %s %s",cl[0],rnu[0],rna[0]);
    ck_assert_str_eq(freesasa_structure_residue_descriptor(s,0),buf);
    for (int i = 0; i < N; ++i) {
        ck_assert_str_eq(freesasa_structure_atom_name(s,i),an[i]);
        ck_assert_str_eq(freesasa_structure_atom_res_name(s,i),rna[i]);
        ck_assert_str_eq(freesasa_structure_atom_res_number(s,i),rnu[i]);
        ck_assert_int_eq(freesasa_structure_atom_chain(s,i),cl[i]);
        sprintf(buf,"%c %s %s %s",cl[i],rnu[i],rna[i],an[i]);
        ck_assert_str_eq(freesasa_structure_atom_descriptor(s,i),buf);
    }
    const freesasa_coord *c = freesasa_structure_xyz(s);
    for (int i = 0; i < N; ++i) {
        const double *xyz = freesasa_coord_i(c, i);
        ck_assert(fabs(xyz[0]+xyz[1]+xyz[2]-3*i) < 1e-10);
    }
    ck_assert(freesasa_structure_n(s) == N);
    ck_assert(freesasa_structure_n_residues(s) == 1);
    
    int first, last;
    ck_assert(freesasa_structure_residue_atoms(s,0,&first,&last) == FREESASA_SUCCESS);
    ck_assert(first == 0 && last == N-1);
}
END_TEST

double a2r(const char *rn, const char *am)
{
    return 1.0;
}

void setup_1ubq(void)
{
    errno = 0;
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(pdb != NULL);
    if (s) freesasa_structure_free(s);
    s = freesasa_structure_from_pdb(pdb,0);
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
    ck_assert(freesasa_structure_n_residues(s) == 76);
    // check at random atom to see that parsing was correct
    ck_assert_str_eq(freesasa_structure_atom_res_name(s,8),"GLN");
    ck_assert_str_eq(freesasa_structure_atom_name(s,8)," N  ");
    ck_assert_str_eq(freesasa_structure_atom_res_number(s,8),"   2");
    ck_assert_int_eq(freesasa_structure_atom_chain(s,8),'A');

    // check coordinates of that random atom
    const freesasa_coord *c = freesasa_structure_xyz(s);
    ck_assert(c != NULL);
    const double *x = freesasa_coord_i(c,8);
    ck_assert(x != NULL);
    ck_assert(fabs(x[0]-26.335+x[1]-27.770+x[2]-3.258) < 1e-10);
}
END_TEST

START_TEST (test_pdb)
{
    const char *file_names[] = {DATADIR "alt_model_twochain.pdb",
                                DATADIR "empty.pdb",
                                DATADIR "empty_model.pdb"};
    const int result_null[] = {0,1,1};
    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (int i = 0; i < 3; ++i) {
        FILE *pdb = fopen(file_names[i],"r");
        ck_assert(pdb != NULL);
        freesasa_structure *s = freesasa_structure_from_pdb(pdb,0);
        if (result_null[i]) ck_assert(s == NULL);
        else ck_assert(s != NULL);
        fclose(pdb);
        freesasa_structure_free(s);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_hydrogen) 
{
    FILE *pdb = fopen(DATADIR "1d3z.pdb","r");
    ck_assert(pdb != NULL);
    freesasa_structure* s = freesasa_structure_from_pdb(pdb,0);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 602);
    freesasa_structure_free(s);
    rewind(pdb);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    s = freesasa_structure_from_pdb(pdb,FREESASA_INCLUDE_HYDROGEN);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 1231);
    freesasa_structure_free(s);
    fclose(pdb);
}
END_TEST

START_TEST (test_hetatm) 
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb","r");
    ck_assert(pdb != NULL);
    freesasa_set_verbosity(FREESASA_V_SILENT); // unknown atoms warnings suppressed
    freesasa_structure* s = freesasa_structure_from_pdb(pdb,FREESASA_INCLUDE_HETATM);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 660);
    freesasa_structure_free(s);
    fclose(pdb);
}
END_TEST

START_TEST (test_structure_array)
{
    FILE *pdb = fopen(DATADIR "1d3z.pdb","r");
    ck_assert(pdb != NULL);
    int n = 0;
    freesasa_structure** ss = freesasa_structure_array(pdb,&n,FREESASA_SEPARATE_MODELS);
        
    ck_assert(ss != NULL);
    ck_assert(n == 10);

    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 602);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    rewind(pdb);
    ss = freesasa_structure_array(pdb,&n,FREESASA_SEPARATE_CHAINS);
    ck_assert(ss != NULL);
    ck_assert(n == 1);
    ck_assert(ss[0] != NULL);
    ck_assert(freesasa_structure_n(ss[0]) == 602);
    free(ss);
    free(ss[0]);
    
    rewind(pdb);
    ss = freesasa_structure_array(pdb,&n,FREESASA_SEPARATE_CHAINS | FREESASA_SEPARATE_MODELS |
                                  FREESASA_INCLUDE_HYDROGEN);
    ck_assert(ss != NULL);
    ck_assert(n == 10);
    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 1231);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    fclose(pdb);
    pdb = fopen(DATADIR "2jo4.pdb", "r");
    ss = freesasa_structure_array(pdb,&n,FREESASA_SEPARATE_MODELS |
                                  FREESASA_INCLUDE_HETATM | FREESASA_INCLUDE_HYDROGEN);
    ck_assert(ss != NULL);
    ck_assert(n == 10);
    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 286*4);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    rewind(pdb);
    pdb = fopen(DATADIR "2jo4.pdb", "r");
    ss = freesasa_structure_array(pdb,&n,FREESASA_SEPARATE_MODELS | FREESASA_SEPARATE_CHAINS |
                                  FREESASA_INCLUDE_HETATM | FREESASA_INCLUDE_HYDROGEN);
    ck_assert(ss != NULL);
    ck_assert(n == 10*4);
    for (int i = 0; i < n; ++i) {
        ck_assert(ss[i] != NULL);
        ck_assert(freesasa_structure_n(ss[i]) == 286);
        freesasa_structure_free(ss[i]);
    }
    free(ss);

    rewind(pdb);
    freesasa_structure *s = freesasa_structure_from_pdb(pdb,FREESASA_INCLUDE_HETATM | 
                                                        FREESASA_INCLUDE_HYDROGEN |
                                                        FREESASA_JOIN_MODELS);
    ck_assert(s != NULL);
    ck_assert(freesasa_structure_n(s) == 286*4*10);
    freesasa_structure_free(s);
    fclose(pdb);
}
END_TEST

Suite* structure_suite() {
    Suite *s = suite_create("Structure");
    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core,setup,teardown);
    tcase_add_test(tc_core, test_structure_api);

    TCase *tc_pdb = tcase_create("PDB");
    tcase_add_test(tc_pdb,test_pdb);
    tcase_add_test(tc_pdb,test_hydrogen);
    tcase_add_test(tc_pdb,test_hetatm);
    tcase_add_test(tc_pdb,test_structure_array);

    TCase *tc_1ubq = tcase_create("1UBQ");
    tcase_add_checked_fixture(tc_1ubq,setup_1ubq,teardown_1ubq);
    tcase_add_test(tc_1ubq,test_structure_1ubq);

    suite_add_tcase(s, tc_core);
    suite_add_tcase(s, tc_pdb);
    suite_add_tcase(s, tc_1ubq);

    return s;
}

