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
#define hyd_O 1.40
#define sulf 2.00
#define n_atom_types 168 //size of array below

const struct atom atoms[] = {
    {"ALA"," C  ",car_C},    {"ALA"," O  ",car_O},    {"ALA"," CA ",ali_C},
    {"ALA"," N  ",ami_N},    {"ALA"," CB ",ali_C}, 
// 5
    {"ARG"," C  ",car_C},    {"ARG"," O  ",car_O},    {"ARG"," CA ",ali_C},
    {"ARG"," N  ",ami_N},    {"ARG"," CB ",ali_C},    {"ARG"," CG ",ali_C},
    {"ARG"," CD ",ali_C},    {"ARG"," CZ ",ali_C},    {"ARG"," NE ",ami_N},
    {"ARG"," NH1",ami_N},    {"ARG"," NH2",ami_N}, 
// 16
    {"ASN"," C  ",car_C},    {"ASN"," O  ",car_O},    {"ASN"," CA ",ali_C},
    {"ASN"," N  ",ami_N},    {"ASN"," CB ",ali_C},    {"ASN"," CG ",ali_C},
    {"ASN"," OD1",car_O},    {"ASN"," ND2",ami_N}, 
// 24
    {"ASP"," C  ",car_C},    {"ASP"," O  ",car_O},    {"ASP"," CA ",ali_C},
    {"ASP"," N  ",ami_N},    {"ASP"," CB ",ali_C},    {"ASP"," CG ",ali_C},
    {"ASP"," OD1",car_O},    {"ASP"," OD2",car_O}, 
// 32
    {"CYS"," C  ",car_C},    {"CYS"," O  ",car_O},    {"CYS"," CA ",ali_C},
    {"CYS"," N  ",ami_N},    {"CYS"," CB ",ali_C},    {"CYS"," SG ",sulf}, 
// 38
    {"GLN"," C  ",car_C},    {"GLN"," O  ",car_O},    {"GLN"," CA ",ali_C},
    {"GLN"," N  ",ami_N},    {"GLN"," CB ",ali_C},    {"GLN"," CG ",ali_C},
    {"GLN"," CD ",ali_C},    {"GLN"," OE1",car_O},    {"GLN"," NE2",ami_N}, 
// 47
    {"GLU"," C  ",car_C},    {"GLU"," O  ",car_O},    {"GLU"," CA ",ali_C},
    {"GLU"," N  ",ami_N},    {"GLU"," CB ",ali_C},    {"GLU"," CG ",ali_C},
    {"GLU"," CD ",ali_C},    {"GLU"," OE1",car_O},    {"GLU"," OE2",car_O}, 
// 56
    {"GLY"," C  ",car_C},    {"GLY"," O  ",car_O},    {"GLY"," CA ",ali_C},
    {"GLY"," N  ",ami_N},    
// 60
    {"HIS"," C  ",car_C},    {"HIS"," O  ",car_O},    {"HIS"," CA ",ali_C},
    {"HIS"," N  ",ami_N},    {"HIS"," CB ",ali_C},    {"HIS"," CG ",aro_C},
    {"HIS"," ND1",ami_N},    {"HIS"," CD2",aro_C},    {"HIS"," NE2",ami_N},
    {"HIS"," CE1",aro_C},
// 70
    {"ILE"," C  ",car_C},    {"ILE"," O  ",car_O},    {"ILE"," CA ",ali_C},
    {"ILE"," N  ",ami_N},    {"ILE"," CB ",ali_C},    {"ILE"," CG1",ali_C},
    {"ILE"," CG2",ali_C},    {"ILE"," CD1",ali_C},
// 78 
    {"LEU"," C  ",car_C},    {"LEU"," O  ",car_O},    {"LEU"," CA ",ali_C},
    {"LEU"," N  ",ami_N},    {"LEU"," CB ",ali_C},    {"LEU"," CG ",ali_C},
    {"LEU"," CD1",ali_C},    {"LEU"," CD2",ali_C},
// 86
    {"LYS"," C  ",car_C},    {"LYS"," O  ",car_O},    {"LYS"," CA ",ali_C},
    {"LYS"," N  ",ami_N},    {"LYS"," CB ",ali_C},    {"LYS"," CG ",ali_C},
    {"LYS"," CG ",ali_C},    {"LYS"," CD ",ali_C},    {"LYS"," CE ",ali_C},
    {"LYS"," NZ ",ami_N},
// 96
    {"MET"," C  ",car_C},    {"MET"," O  ",car_O},    {"MET"," CA ",ali_C},
    {"MET"," N  ",ami_N},    {"MET"," CB ",ali_C},    {"MET"," CG ",ali_C},
    {"MET"," SD ",sulf},     {"MET"," CE ",ali_C},
// 104
    {"PHE"," C  ",car_C},    {"PHE"," O  ",car_O},    {"PHE"," CA ",ali_C},
    {"PHE"," N  ",ami_N},    {"PHE"," CB ",ali_C},    {"PHE"," CG ",aro_C},
    {"PHE"," CD1",aro_C},    {"PHE"," CD2",aro_C},    {"PHE"," CE1",aro_C},
    {"PHE"," CE2",aro_C},    {"PHE"," CZ ",aro_C},
// 115
    {"PRO"," C  ",car_C},    {"PRO"," O  ",car_O},    {"PRO"," CA ",ali_C},
    {"PRO"," N  ",ami_N},    {"PRO"," CB ",aro_C},    {"PRO"," CG ",aro_C},
    {"PRO"," CD ",aro_C},
// 122
    {"SER"," C  ",car_C},    {"SER"," O  ",car_O},    {"SER"," CA ",ali_C},
    {"SER"," N  ",ami_N},    {"SER"," CB ",ali_C},    {"SER"," OG ",hyd_O},
// 128
    {"THR"," C  ",car_C},    {"THR"," O  ",car_O},    {"THR"," CA ",ali_C},
    {"THR"," N  ",ami_N},    {"THR"," CB ",ali_C},    {"THR"," OG1",hyd_O},
    {"THR"," CG2",ali_C},
// 135
    {"TRP"," C  ",car_C},    {"TRP"," O  ",car_O},    {"TRP"," CA ",ali_C},
    {"TRP"," N  ",ami_N},    {"TRP"," CB ",ali_C},    {"TRP"," CG ",aro_C},
    {"TRP"," CD1",aro_C},    {"TRP"," CD2",aro_C},    {"TRP"," NE1",ami_N},
    {"TRP"," CE2",aro_C},    {"TRP"," CE3",aro_C},    {"TRP"," CZ2",aro_C},
    {"TRP"," CZ3",aro_C},    {"TRP"," CH2",aro_C},
// 149
    {"TYR"," C  ",car_C},    {"TYR"," O  ",car_O},    {"TYR"," CA ",ali_C},
    {"TYR"," N  ",ami_N},    {"TYR"," CB ",ali_C},    {"TYR"," CG ",aro_C},
    {"TYR"," CD1",aro_C},    {"TYR"," CD2",aro_C},    {"TYR"," CE1",aro_C},
    {"TYR"," CE2",aro_C},    {"TYR"," CZ ",aro_C},    {"TYR"," OH ",hyd_O},
// 161
    {"VAL"," C  ",car_C},    {"VAL"," O  ",car_O},    {"VAL"," CA ",ali_C},
    {"VAL"," N  ",ami_N},    {"VAL"," CB ",ali_C},    {"VAL"," CG1",ali_C},
    {"VAL"," CG2",ali_C},
// 168
};

START_TEST (test_radius)
{
    char buf[50];
    struct atom a;
    for (int i = 0; i < n_atom_types; ++i) {
        a = atoms[i];
        double ret = sasalib_classify_radius(a.a,a.b);
        sprintf(buf,"%s %s ret=%f ref=%f",a.a,a.b,ret,a.radius);
        ck_assert_msg(fabs(ret - a.radius) < 1e-10,buf);
    }
    //printf("%s %s\n",a.a,a.b);
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
