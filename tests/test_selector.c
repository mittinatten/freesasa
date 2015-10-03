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

#include <string.h>
#include <math.h>
#include <stdio.h>
#include <freesasa.h>
#include <check.h>

#define N 8
freesasa_structure *structure;
freesasa_result *result;
freesasa_strvp *svp;
double radii[N] = {1.0,1.2,1.4,1.6,1.8,2,2.2,2.4};
const char *name[N] = { " CA ", " O  ", " N  ", " SD ", " CB ", " OXT", " SE ", " P  "};
const char *resn[N] = {  "ALA",  "ALA",  "ARG",  "MET",  "VAL",  "GLU",  "SEC",  " DC"};
const char *resi[N] = { "   1", "   1", "   2", "   3", "   4", "   1", "   2", "   4"};
const char chain[N] = {    'A',    'A',    'A',    'A',    'A',    'B',    'B',    'B'};
const int symb_O[N] = {      0,      1,      0,      0,      0,     1,       0,      0};
const int symb_C[N] = {      1,      0,      0,      0,      1,     0,       0,      0};
const int symb_N[N] = {      0,      0,      1,      0,      0,     0,       0,      0};
const int symb_S[N] = {      0,      0,      0,      1,      0,     0,       0,      0};
const int symb_SE[N]= {      0,      0,      0,      0,      0,     0,       1,      0};
const int symb_P[N] = {      0,      0,      0,      0,      0,     0,       0,      1};
const int resn_A[N] = {      1,      1,      0,      0,      0,     0,       0,      0};
const int resn_R[N] = {      0,      0,      1,      0,      0,     0,       0,      0};
const int resn_M[N] = {      0,      0,      0,      1,      0,     0,       0,      0};
const int resn_V[N] = {      0,      0,      0,      0,      1,     0,       0,      0};
const int resn_Q[N] = {      0,      0,      0,      0,      0,     1,       0,      0};
const int resn_U[N] = {      0,      0,      0,      0,      0,     0,       1,      0};
const int resn_DC[N]= {      0,      0,      0,      0,      0,     0,       0,      1};
const int name_CA[N]= {      1,      0,      0,      0,      0,     0,       0,      0};
const int name_O[N] = {      0,      1,      0,      0,      0,     0,       0,      0};
const int name_OXT[N]={      0,      0,      0,      0,      0,     1,       0,      0};

int float_eq(double a, double b, double tolerance) {
    if (fabs(a-b) < tolerance) return 1;
    printf("floats not equal: a = %f, b = %f, diff = %f, tolerance = %f\n",
           a, b, fabs(a-b), tolerance);
    fflush(stdout);
    return 0;
}

static void setup(void) 
{
    structure = freesasa_structure_new();
    for (int i = 0; i < N; ++i) {
        freesasa_structure_add_atom(structure,name[i],resn[i],resi[i],chain[i],i*10,0,0);
    }
    result = freesasa_calc_structure(structure, radii, NULL);
}


static void teardown(void) 
{
    freesasa_structure_free(structure);
    freesasa_result_free(result);
    freesasa_strvp_free(svp);
}

static double addup(const int *sel, const freesasa_result *res)
{
    double sum = 0;
    for (int i = 0; i < res->n_atoms; ++i) 
        sum += sel[i]*res->sasa[i];
    return sum;
}

START_TEST (test_name)
{
    const char *commands[] = {"c1, name ca+o",
                              "c2, name ca",
                              "c3, name oxt",
                              "c4, name ca AND name o",
                              "c5, name ca OR  name o"};
    svp = freesasa_select_area(commands,5,structure,result);
    ck_assert_ptr_ne(svp,NULL);
    ck_assert_ptr_ne(svp->value,NULL);
    ck_assert_ptr_ne(svp->string,NULL);
    ck_assert_int_eq(svp->n,5);
    ck_assert(svp->value[0] > 5); // check that it's non-zero
    ck_assert(float_eq(svp->value[0], addup(name_CA,result) + addup(name_O,result), 1e-10));
    ck_assert(float_eq(svp->value[0], svp->value[4], 1e-10));
    ck_assert(float_eq(svp->value[1], addup(name_CA,result), 1e-10));
    ck_assert(float_eq(svp->value[2], addup(name_OXT,result), 1e-10));
    ck_assert(float_eq(svp->value[3], 0, 1e-10));
}
END_TEST

START_TEST (test_symbol)
{
    const char *commands[] = {"c1, symbol o+c",
                              "c2, symbol O",
                              "c3, symbol C",
                              "c4, symbol O AND symbol C",
                              "c5, symbol O OR symbol C"};
    svp = freesasa_select_area(commands,5,structure,result);
    ck_assert_ptr_ne(svp,NULL);
    ck_assert_ptr_ne(svp->value,NULL);
    ck_assert_ptr_ne(svp->string,NULL);
    ck_assert_int_eq(svp->n,5);
    ck_assert(svp->value[0] > 5); //just to check that it's non-zero
    ck_assert(float_eq(svp->value[0], addup(symb_O,result) + addup(symb_C,result), 1e-10));
    ck_assert(float_eq(svp->value[0], svp->value[4], 1e-10));
    ck_assert(float_eq(svp->value[1], addup(symb_O,result), 1e-10));
    ck_assert(float_eq(svp->value[2], addup(symb_C,result), 1e-10));
    ck_assert(float_eq(svp->value[3], 0, 1e-10));
    ck_assert_ptr_ne(svp->string[0], NULL);
    ck_assert_str_eq(svp->string[0], "c1");
    ck_assert_str_eq(svp->string[1], "c2");
    ck_assert_str_eq(svp->string[2], "c3");
}
END_TEST

START_TEST (test_resn)
{
    const char *commands[] = {"c1, resn ala+arg",
                              "c2, resn ala",
                              "c3, resn arg",
                              "c4, resn ala AND resn arg",
                              "c5, resn ala OR  resn arg"};
    svp = freesasa_select_area(commands,5,structure,result);
    ck_assert_ptr_ne(svp,NULL);
    ck_assert_ptr_ne(svp->value,NULL);
    ck_assert_ptr_ne(svp->string,NULL);
    ck_assert_int_eq(svp->n,5);
    ck_assert(svp->value[0] > 5);
    ck_assert(float_eq(svp->value[0], addup(resn_A,result) + addup(resn_R,result), 1e-10));
    ck_assert(float_eq(svp->value[0], svp->value[4], 1e-10));
    ck_assert(float_eq(svp->value[1], addup(resn_A,result), 1e-10));
    ck_assert(float_eq(svp->value[2], addup(resn_R,result), 1e-10));
    ck_assert(float_eq(svp->value[3], 0, 1e-10));
}
END_TEST

Suite *selector_suite() {
    Suite *s = suite_create("Selector");

    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core,setup,teardown);
    tcase_add_test(tc_core, test_name);
    tcase_add_test(tc_core, test_symbol);
    tcase_add_test(tc_core, test_resn);

    suite_add_tcase(s, tc_core);

    return s;
}
