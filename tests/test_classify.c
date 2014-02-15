#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <errno.h>
#include <check.h>
#include <classify.h>

struct atom {
    const char *a;
    const char *b;
    const double radius;
};
#define ali_C 2.00
#define aro_C 1.75
#define car_O 1.40
#define car_C 1.55
#define ami_N 1.55
#define n_atom_types 31 //size of array below

const struct atom atoms[] = {
    {"ALA"," C  ",car_C},    {"ALA"," O  ",car_O},    {"ALA"," CA ",ali_C},
    {"ALA"," N  ",ami_N},    {"ALA"," CB ",ali_C},
    {"ARG"," C  ",car_C},    {"ARG"," O  ",car_O},    {"ARG"," CA ",ali_C},
    {"ARG"," N  ",ami_N},    {"ARG"," CB ",ali_C},    {"ARG"," CG ",ali_C},
    {"ARG"," CD ",ali_C},    {"ARG"," CZ ",ali_C},    {"ARG"," NE ",ami_N},
    {"ARG"," NH1",ami_N},    {"ARG"," NH2",ami_N},
    {"ASN"," C  ",car_C},    {"ASN"," O  ",car_O},    {"ASN"," CA ",ali_C},
    {"ASN"," N  ",ami_N},    {"ASN"," CB ",ali_C},    {"ASN"," CG ",ali_C},
    {"ASN"," OD1",car_O},    {"ASN"," ND2",ami_N},
    {"ASP"," C  ",car_C},    {"ASP"," O  ",car_O},    {"ASP"," CA ",ali_C},
    {"ASP"," N  ",ami_N},    {"ASP"," CB ",ali_C},    {"ASP"," CG ",ali_C},
    {"ASN"," OD1",car_O},    {"ASN"," OD2",car_O},
};

START_TEST (test_radius)
{
    char buf[50];
    for (int i = 0; i < n_atom_types; ++i) {
        const struct atom a = atoms[i];
        double ret = sasalib_classify_radius(a.a,a.b);
        sprintf(buf,"%s %s ret=%f ref=%f",a.a,a.b,ret,a.radius);
        ck_assert_msg(fabs(ret - a.radius) < 1e-10,buf);
    }
}
END_TEST

START_TEST (test_residue) 
{
    int nrt = sasalib_classify_nresiduetypes();
    const char **res = (const char**) malloc(nrt*sizeof(char*));
    for (int i = 0; i < nrt; ++i) {
	res[i] = sasalib_classify_residue2str(i);
	ck_assert_int_eq(sasalib_classify_residue(res[i]),i);
    }
    ck_assert(sasalib_classify_residue("ALA") < 20);
    ck_assert(sasalib_classify_residue("ARG") < 20);
    ck_assert(sasalib_classify_residue("ASN") < 20);
    ck_assert(sasalib_classify_residue("ASP") < 20);
    ck_assert(sasalib_classify_residue("CYS") < 20);
    ck_assert(sasalib_classify_residue("GLN") < 20);
    ck_assert(sasalib_classify_residue("GLU") < 20);
    ck_assert(sasalib_classify_residue("GLY") < 20);
    ck_assert(sasalib_classify_residue("HIS") < 20);
    ck_assert(sasalib_classify_residue("ILE") < 20);
    ck_assert(sasalib_classify_residue("LEU") < 20);
    ck_assert(sasalib_classify_residue("LYS") < 20);
    ck_assert(sasalib_classify_residue("MET") < 20);
    ck_assert(sasalib_classify_residue("PHE") < 20);
    ck_assert(sasalib_classify_residue("PRO") < 20);
    ck_assert(sasalib_classify_residue("SER") < 20);
    ck_assert(sasalib_classify_residue("THR") < 20);
    ck_assert(sasalib_classify_residue("TRP") < 20);
    ck_assert(sasalib_classify_residue("TYR") < 20);
    ck_assert(sasalib_classify_residue("VAL") < 20);
}
END_TEST

Suite* classify_suite() 
{
    Suite *s = suite_create("Classify");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core,test_radius);
    tcase_add_test(tc_core,test_residue);

    
    suite_add_tcase(s,tc_core);

    return s;
}
