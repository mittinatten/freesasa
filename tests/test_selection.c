#include <string.h>
#include <math.h>
#include <stdio.h>
#include <freesasa.h>
#include <check.h>
#include "tools.h"

#define N 8

freesasa_structure *structure;
freesasa_result *result;
freesasa_strvp *svp;
double radii[N] = {1.0,1.2,1.4,1.6,1.8,2,2.2,2.4};
const char *name[N] = { " CA ", " O  ", " N  ", " SD ", " CB ", " OXT", "SE  ", " P  "};
const char *resn[N] = {  "ALA",  "ALA",  "ARG",  "MET",  "VAL",  "GLU",  "SEC",  " DC"};
const char *resi[N] = { "   1", "   1", "   2", "   3", "   4", "   1", "   2", "   4"};
const char chain[N] = {    'A',    'A',    'A',    'A',    'A',    'B',    'B',    'B'};
const int symb_O[N] = {      0,      1,      0,      0,      0,      1,      0,      0};
const int symb_C[N] = {      1,      0,      0,      0,      1,      0,      0,      0};
const int symb_N[N] = {      0,      0,      1,      0,      0,      0,      0,      0};
const int symb_S[N] = {      0,      0,      0,      1,      0,      0,      0,      0};
const int symb_SE[N]= {      0,      0,      0,      0,      0,      0,      1,      0};
const int symb_P[N] = {      0,      0,      0,      0,      0,      0,      0,      1};
const int resn_A[N] = {      1,      1,      0,      0,      0,      0,      0,      0};
const int resn_R[N] = {      0,      0,      1,      0,      0,      0,      0,      0};
const int resn_M[N] = {      0,      0,      0,      1,      0,      0,      0,      0};
const int resn_V[N] = {      0,      0,      0,      0,      1,      0,      0,      0};
const int resn_Q[N] = {      0,      0,      0,      0,      0,      1,      0,      0};
const int resn_U[N] = {      0,      0,      0,      0,      0,      0,      1,      0};
const int resn_DC[N]= {      0,      0,      0,      0,      0,      0,      0,      1};
const int name_CA[N]= {      1,      0,      0,      0,      0,      0,      0,      0};
const int name_O[N] = {      0,      1,      0,      0,      0,      0,      0,      0};
const int name_OXT[N]={      0,      0,      0,      0,      0,      1,      0,      0};
const int resi_1r4[N]={      1,      1,      1,      1,      1,      1,      1,      1};
const int resi_1[N] = {      1,      1,      0,      0,      0,      1,      0,      0};
const int resi_2r4[N]={      0,      0,      1,      1,      1,      0,      1,      1};
const int chain_A[N]= {      1,      1,      1,      1,      1,      0,      0,      0};
const int chain_B[N]= {      0,      0,      0,      0,      0,      1,      1,      1};
char selection_name[100][FREESASA_MAX_SELECTION_NAME+1];
double value[100];

static void
test_select(const char **command,int n_commands) 
{
    for (int i = 0; i < n_commands; ++i) 
        ck_assert_int_eq(freesasa_select_area(command[i],selection_name[i],value+i,structure,result),
                         FREESASA_SUCCESS);
}

static void setup(void) 
{
    structure = freesasa_structure_new();
    for (int i = 0; i < N; ++i) {
        freesasa_structure_add_atom(structure,name[i],resn[i],resi[i],chain[i],i*10,0,0);
    }
    result = freesasa_calc_structure(structure, NULL);
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
                              "c5, name ca OR  name o",
                              "c6, name ca+o+oxt"};
    test_select(commands,6);
    ck_assert(value[0] > 5); // check that it's non-zero
    ck_assert(float_eq(value[0], addup(name_CA,result) + addup(name_O,result), 1e-10));
    ck_assert(float_eq(value[0], value[4], 1e-10));
    ck_assert(float_eq(value[1], addup(name_CA,result), 1e-10));
    ck_assert(float_eq(value[2], addup(name_OXT,result), 1e-10));
    ck_assert(float_eq(value[3], 0, 1e-10));
    ck_assert(float_eq(value[5], 
                       addup(name_CA,result) + addup(name_O,result) + addup(name_OXT,result),
                       1e-10));

    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_select_area("c1, name ABCDE",selection_name[0],value,structure,result),
                     FREESASA_WARN);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_symbol)
{
    const char *commands[] = {"c1, symbol o+c",
                              "c2, symbol O",
                              "c3, symbol C",
                              "c4, symbol O AND symbol C",
                              "c5, symbol O OR symbol C",
                              "c6, symbol O+C+SE",
                              "c7, symbol SE",
                              "c8, symbol O+C+SE and not symbol se"};
    test_select(commands,8);
    ck_assert_str_eq(selection_name[0], "c1");
    ck_assert_str_eq(selection_name[1], "c2");
    ck_assert_str_eq(selection_name[2], "c3");
    ck_assert(value[0] > 5); //just to check that it's non-zero
    ck_assert(float_eq(value[0], addup(symb_O,result) + addup(symb_C,result), 1e-10));
    ck_assert(float_eq(value[1], addup(symb_O,result), 1e-10));
    ck_assert(float_eq(value[2], addup(symb_C,result), 1e-10));
    ck_assert(float_eq(value[3], 0, 1e-10));
    ck_assert(float_eq(value[4], value[0], 1e-10));
    ck_assert(float_eq(value[5], 
                       addup(symb_O,result) + addup(symb_C,result) + addup(symb_SE,result),
                       1e-10));
    ck_assert(float_eq(value[6], addup(symb_SE,result), 1e-10));
    ck_assert(float_eq(value[7], value[0], 1e-10));

    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_select_area("c1, symbol ABC",selection_name[0],value,structure,result),
                     FREESASA_WARN);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_resn)
{
    const char *commands[] = {"c1, resn ala+arg",
                              "c2, resn ala",
                              "c3, resn arg",
                              "c4, resn ala AND resn arg",
                              "c5, resn ala OR  resn arg",
                              "c6, resn ala+arg AND NOT resn arg"};
    test_select(commands,6);
    ck_assert(value[0] > 5);
    ck_assert(float_eq(value[0], addup(resn_A,result) + addup(resn_R,result), 1e-10));
    ck_assert(float_eq(value[1], addup(resn_A,result), 1e-10));
    ck_assert(float_eq(value[2], addup(resn_R,result), 1e-10));
    ck_assert(float_eq(value[3], 0, 1e-10));
    ck_assert(float_eq(value[4], value[0], 1e-10));
    ck_assert(float_eq(value[5], value[1], 1e-10));

    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_select_area("c1, resn ABCD",selection_name[0],value,structure,result),
                     FREESASA_WARN);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_resi)
{
    const char *commands[] = {"c1, resi 1+2-4",
                              "c2, resi 2-4",
                              "c3, resi 1",
                              "c4, resi 1 AND resi 2-4",
                              "c5, resi 1 OR  resi 2-4",
                              "c6, resi 1-2+2-4",
                              "c7, resi 1+2-4+3",
                              "c8, resi 1-2+7+9+3-5+100",
                              "c9, resi 1-4 AND NOT resi 2-4"};
    freesasa_set_verbosity(FREESASA_V_SILENT);
    test_select(commands,9);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    ck_assert(value[0] > 5);
    ck_assert(float_eq(value[0], addup(resi_1,result) + addup(resi_2r4,result), 1e-10));
    ck_assert(float_eq(value[1], addup(resi_2r4,result), 1e-10));
    ck_assert(float_eq(value[2], addup(resi_1,result), 1e-10));
    ck_assert(float_eq(value[3], 0, 1e-10));
    ck_assert(float_eq(value[4], value[0], 1e-10));
    ck_assert(float_eq(value[5], value[0], 1e-10));
    ck_assert(float_eq(value[6], value[0], 1e-10));
    ck_assert(float_eq(value[7], value[0], 1e-10));
    ck_assert(float_eq(value[8], value[2], 1e-10));

    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_select_area("c1, resi A",selection_name[0],value,structure,result),
                     FREESASA_WARN);
    ck_assert_int_eq(freesasa_select_area("c1, resi 1A",selection_name[0],value,structure,result),
                     FREESASA_FAIL); // This is a syntax error
    ck_assert_int_eq(freesasa_select_area("c1, resi 1-A",selection_name[0],value,structure,result),
                     FREESASA_WARN);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_chain)
{
    const char *commands[] = {"c1, chain A+B",
                              "c2, chain A",
                              "c3, chain B",
                              "c4, chain A AND chain B",
                              "c5, chain A OR chain B",
                              "c6, chain A-B",
                              "c7, chain A-B AND NOT chain A"};
    test_select(commands,7);
    ck_assert(value[0] > 5);
    ck_assert(float_eq(value[0], addup(chain_A,result) + addup(chain_B,result), 1e-10));
    ck_assert(float_eq(value[0], value[4], 1e-10));
    ck_assert(float_eq(value[0], value[5], 1e-10));
    ck_assert(float_eq(value[1], addup(chain_A,result), 1e-10));
    ck_assert(float_eq(value[2], addup(chain_B,result), 1e-10));
    ck_assert(float_eq(value[3], 0, 1e-10));
    ck_assert(float_eq(value[6], addup(chain_B,result), 1e-10));

    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(freesasa_select_area("c1, chain AA",selection_name[0],value,structure,result),
                     FREESASA_WARN);
    ck_assert_int_eq(freesasa_select_area("c1, chain &",selection_name[0],value,structure,result),
                     FREESASA_FAIL);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST (test_syntax_error)
{
    double a;
    int i;
    char s[FREESASA_MAX_SELECTION_NAME+1];
    char s2[2*FREESASA_MAX_SELECTION_NAME+1], c='a';
    const char *err[] = {"","a","a,","a,b",
                         // no selection arg
                         "a,resi","a,resn","a,name","a,symbol","a,chain",
                         // no name of selection
                         ",resn ala",",resi 1",",name ca",", symbol c",",chain a",
                         "resn ala","resi 1","name ca","symbol c","chain a",
                         // comma wrong place, no name
                         "resn, ala","resi, 1","name, ca","symbol, c","chain, a",
                         // ranges (-) used where not allowed
                         "a, resn ala-arg", "a, name ca-cb","a, symbol c-o",
                         "a, resi 1-2-3",
                         // trailing +-
                         "a, resn ala+", "a, resn ala+arg+", "a, resi 1-",
                         "a, resi 1-", "a, resi 1-2+","a, resi 1+2-5+",
                         // boolean operators
                         "a, (resn ala) AND","a,(resn ala) OR","a,(resn ala) OR NOT",
                         "a, (resn ala) AND arg","a,(resn ala) OR arg",
                         "a, (resn ala) OR NOT arg",
                         "a, ala OR resn arg",
    };
    int n = sizeof(err)/sizeof(char*);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (i = 0; i < n; ++i) {
        ck_assert_int_eq(freesasa_select_area(err[i],s,&a,structure,result), FREESASA_FAIL);
    }

    // check that really long strings are truncated properly
    for (i = 0; i < FREESASA_MAX_SELECTION_NAME+5; ++i) {
        s2[i] = c+(i%25);
    }
    strcpy(s2+i, ", resn ala");
    ck_assert_int_eq(freesasa_select_area(s2, s, &a, structure, result), FREESASA_SUCCESS);
    ck_assert_int_eq(strlen(s), FREESASA_MAX_SELECTION_NAME);
    ck_assert(strncmp(s2, s, strlen(s)) == 0);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
} END_TEST

  // check that some complex constructs don't cause errors, but no checks for valid selection
START_TEST (test_complex_syntax)
{
    double a;
    char s[FREESASA_MAX_SELECTION_NAME+1];
    const char *stmt[] = {"a, (resn ala AND resi 1-3) OR (NOT chain A+B AND (symbol C OR symbol O))",
                          "a, NOT symbol SE+C AND NOT resi 5-7+1+6-8+100+200+10-11"
    };
    int n = sizeof(stmt)/sizeof(char*);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (int i = 0; i < n; ++i) {
        ck_assert_int_eq(freesasa_select_area(stmt[i],s,&a,structure,result),FREESASA_SUCCESS);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
} END_TEST

  // This test fails because the Bison-generated parser breaks
  // when memory allocation fails
/* 
START_TEST (test_memerr) 
{
    set_fail_after(0);
    freesasa_parameters p = freesasa_default_parameters;
    p.shrake_rupley_n_points = 10;
    FILE *file = fopen(DATADIR "1ubq.pdb", "r");
    freesasa_structure *s = freesasa_structure_from_pdb(file, NULL, 0);
    freesasa_result *result = freesasa_calc_structure(s, &p);
    double area;
    char name[FREESASA_MAX_SELECTION_NAME];
    fclose(file);

    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (int i = 1; i < 17; ++i) { 
    set_fail_after(i);
        // this expression should come across most allocations
        int ret = freesasa_select_area("s, resn ALA and not resi 1-20", name, &area, s, result);
        set_fail_after(0);
        ck_assert_int_eq(ret, FREESASA_FAIL);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
    freesasa_result_free(result);
    freesasa_structure_free(s);
}
END_TEST
*/

Suite *selector_suite() {
    Suite *s = suite_create("Selector");

    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core,setup,teardown);
    tcase_add_test(tc_core, test_name);
    tcase_add_test(tc_core, test_symbol);
    tcase_add_test(tc_core, test_resn);
    tcase_add_test(tc_core, test_resi);
    tcase_add_test(tc_core, test_chain);
    //tcase_add_test(tc_core, test_memerr);
    
    TCase *tc_syntax = tcase_create("Syntax");
    // just to avoid passing NULL pointers
    tcase_add_checked_fixture(tc_syntax,setup,teardown);
    tcase_add_test(tc_syntax, test_syntax_error);
    tcase_add_test(tc_syntax, test_complex_syntax);

    suite_add_tcase(s, tc_core);
    suite_add_tcase(s, tc_syntax);

    return s;
}
