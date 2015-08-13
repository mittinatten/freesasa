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
#include <classify.h>

struct atom {
    const char *a;
    const char *b;
    const double radius;
    const int class;
};
#define ali_C 2.00
#define aro_C 1.75
#define car_O 1.40
#define car_C 1.55
#define ami_N 1.55
#define hyd_O 1.40
#define sulf 2.00
#define sele 1.90
#define unk_P 1.5
#define r_unk 0.0
#define APO FREESASA_APOLAR
#define POL FREESASA_POLAR
#define NUC FREESASA_NUCLEICACID
#define UNK FREESASA_CLASS_UNKNOWN

#define naa 20
const char *aa[naa] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                       "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                       "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
#define n_other 3
const char *other_aa[n_other] = {"GLX", "ASX", "CSE"};

#define n_nuc 13
const char *nuc[n_nuc] = {"DA","DC","DG","DT","DU","DI","A","C","G","T","U","I","N"};


#define n_atom_types 218
const struct atom atoms[n_atom_types] = {
    {"ALA"," C  ",car_C,POL}, {"ALA"," O  ",car_O,POL}, {"ALA"," CA ",ali_C,APO},
    {"ALA"," N  ",ami_N,POL}, {"ALA"," CB ",ali_C,APO}, {"ALA"," X  ",r_unk,UNK},
// 6
    {"ARG"," C  ",car_C,POL}, {"ARG"," O  ",car_O,POL}, {"ARG"," CA ",ali_C,APO},
    {"ARG"," N  ",ami_N,POL}, {"ARG"," CB ",ali_C,APO}, {"ARG"," CG ",ali_C,APO},
    {"ARG"," CD ",ali_C,APO}, {"ARG"," CZ ",ali_C,APO}, {"ARG"," NE ",ami_N,POL},
    {"ARG"," NH1",ami_N,POL}, {"ARG"," NH2",ami_N,POL}, {"ARG"," X  ",r_unk,UNK},
// 18
    {"ASN"," C  ",car_C,POL}, {"ASN"," O  ",car_O,POL}, {"ASN"," CA ",ali_C,APO},
    {"ASN"," N  ",ami_N,POL}, {"ASN"," CB ",ali_C,APO}, {"ASN"," CG ",car_C,POL},
    {"ASN"," OD1",car_O,POL}, {"ASN"," ND2",ami_N,POL}, {"ASN"," Y  ",r_unk,UNK},
// 27
    {"ASP"," C  ",car_C,POL}, {"ASP"," O  ",car_O,POL}, {"ASP"," CA ",ali_C,APO},
    {"ASP"," N  ",ami_N,POL}, {"ASP"," CB ",ali_C,APO}, {"ASP"," CG ",car_C,POL},
    {"ASP"," OD1",car_O,POL}, {"ASP"," OD2",car_O,POL}, {"ASP"," Y  ",r_unk,UNK},
// 36
    {"CYS"," C  ",car_C,POL}, {"CYS"," O  ",car_O,POL}, {"CYS"," CA ",ali_C,APO},
    {"CYS"," N  ",ami_N,POL}, {"CYS"," CB ",ali_C,APO}, {"CYS"," SG ",sulf,POL},
    {"CYS"," X  ",r_unk,UNK},
// 43
    {"GLN"," C  ",car_C,POL}, {"GLN"," O  ",car_O,POL}, {"GLN"," CA ",ali_C,APO},
    {"GLN"," N  ",ami_N,POL}, {"GLN"," CB ",ali_C,APO}, {"GLN"," CG ",ali_C,APO},
    {"GLN"," CD ",car_C,POL}, {"GLN"," OE1",car_O,POL}, {"GLN"," NE2",ami_N,POL},
    {"GLN"," Y  ",r_unk,UNK},
// 53
    {"GLU"," C  ",car_C,POL}, {"GLU"," O  ",car_O,POL}, {"GLU"," CA ",ali_C,APO},
    {"GLU"," N  ",ami_N,POL}, {"GLU"," CB ",ali_C,APO}, {"GLU"," CG ",ali_C,APO},
    {"GLU"," CD ",car_C,POL}, {"GLU"," OE1",car_O,POL}, {"GLU"," OE2",car_O,POL},
    {"GLU"," Y  ",r_unk,UNK},
// 63
    {"GLY"," C  ",car_C,POL}, {"GLY"," O  ",car_O,POL}, {"GLY"," CA ",ali_C,APO},
    {"GLY"," N  ",ami_N,POL}, {"GLY"," X  ",r_unk,UNK},
// 68
    {"HIS"," C  ",car_C,POL}, {"HIS"," O  ",car_O,POL}, {"HIS"," CA ",ali_C,APO},
    {"HIS"," N  ",ami_N,POL}, {"HIS"," CB ",ali_C,APO}, {"HIS"," CG ",aro_C,APO},
    {"HIS"," ND1",ami_N,POL}, {"HIS"," CD2",aro_C,APO}, {"HIS"," NE2",ami_N,POL},
    {"HIS"," CE1",aro_C,APO}, {"HIS"," X  ",r_unk,UNK},
// 79
    {"ILE"," C  ",car_C,POL}, {"ILE"," O  ",car_O,POL}, {"ILE"," CA ",ali_C,APO},
    {"ILE"," N  ",ami_N,POL}, {"ILE"," CB ",ali_C,APO}, {"ILE"," CG1",ali_C,APO},
    {"ILE"," CG2",ali_C,APO}, {"ILE"," CD1",ali_C,APO}, {"ILE"," X  ",r_unk,UNK},
// 88
    {"LEU"," C  ",car_C,POL}, {"LEU"," O  ",car_O,POL}, {"LEU"," CA ",ali_C,APO},
    {"LEU"," N  ",ami_N,POL}, {"LEU"," CB ",ali_C,APO}, {"LEU"," CG ",ali_C,APO},
    {"LEU"," CD1",ali_C,APO}, {"LEU"," CD2",ali_C,APO}, {"LEU"," X  ",r_unk,UNK},
// 97
    {"LYS"," C  ",car_C,POL}, {"LYS"," O  ",car_O,POL}, {"LYS"," CA ",ali_C,APO},
    {"LYS"," N  ",ami_N,POL}, {"LYS"," CB ",ali_C,APO}, {"LYS"," CG ",ali_C,APO},
    {"LYS"," CG ",ali_C,APO}, {"LYS"," CD ",ali_C,APO}, {"LYS"," CE ",ali_C,APO},
    {"LYS"," NZ ",ami_N,POL}, {"LYS"," X  ",r_unk,UNK},
// 108
    {"MET"," C  ",car_C,POL}, {"MET"," O  ",car_O,POL}, {"MET"," CA ",ali_C,APO},
    {"MET"," N  ",ami_N,POL}, {"MET"," CB ",ali_C,APO}, {"MET"," CG ",ali_C,APO},
    {"MET"," SD ",sulf},  {"MET"," CE ",ali_C,APO}, {"MET"," X  ",r_unk,UNK},
// 117
    {"PHE"," C  ",car_C,POL}, {"PHE"," O  ",car_O,POL}, {"PHE"," CA ",ali_C,APO},
    {"PHE"," N  ",ami_N,POL}, {"PHE"," CB ",ali_C,APO}, {"PHE"," CG ",aro_C,APO},
    {"PHE"," CD1",aro_C,APO}, {"PHE"," CD2",aro_C,APO}, {"PHE"," CE1",aro_C,APO},
    {"PHE"," CE2",aro_C,APO}, {"PHE"," CZ ",aro_C,APO}, {"PHE"," X  ",r_unk,UNK},
// 129
    {"PRO"," C  ",car_C,POL}, {"PRO"," O  ",car_O,POL}, {"PRO"," CA ",ali_C,APO},
    {"PRO"," N  ",ami_N,POL}, {"PRO"," CB ",aro_C,APO}, {"PRO"," CG ",aro_C,APO},
    {"PRO"," CD ",aro_C,APO}, {"PRO"," X  ",r_unk,UNK},
// 137
    {"SER"," C  ",car_C,POL}, {"SER"," O  ",car_O,POL}, {"SER"," CA ",ali_C,APO},
    {"SER"," N  ",ami_N,POL}, {"SER"," CB ",ali_C,APO}, {"SER"," OG ",hyd_O,POL},
    {"SER"," X  ",r_unk,UNK},
// 144
    {"THR"," C  ",car_C,POL}, {"THR"," O  ",car_O,POL}, {"THR"," CA ",ali_C,APO},
    {"THR"," N  ",ami_N,POL}, {"THR"," CB ",ali_C,APO}, {"THR"," OG1",hyd_O,POL},
    {"THR"," CG2",ali_C,APO}, {"THR"," X  ",r_unk,UNK},
// 152
    {"TRP"," C  ",car_C,POL}, {"TRP"," O  ",car_O,POL}, {"TRP"," CA ",ali_C,APO},
    {"TRP"," N  ",ami_N,POL}, {"TRP"," CB ",ali_C,APO}, {"TRP"," CG ",aro_C,APO},
    {"TRP"," CD1",aro_C,APO}, {"TRP"," CD2",aro_C,APO}, {"TRP"," NE1",ami_N,POL},
    {"TRP"," CE2",aro_C,APO}, {"TRP"," CE3",aro_C,APO}, {"TRP"," CZ2",aro_C,APO},
    {"TRP"," CZ3",aro_C,APO}, {"TRP"," CH2",aro_C,APO}, {"TRP"," X  ",r_unk,UNK},
// 167
    {"TYR"," C  ",car_C,POL}, {"TYR"," O  ",car_O,POL}, {"TYR"," CA ",ali_C,APO},
    {"TYR"," N  ",ami_N,POL}, {"TYR"," CB ",ali_C,APO}, {"TYR"," CG ",aro_C,APO},
    {"TYR"," CD1",aro_C,APO}, {"TYR"," CD2",aro_C,APO}, {"TYR"," CE1",aro_C,APO},
    {"TYR"," CE2",aro_C,APO}, {"TYR"," CZ ",aro_C,APO}, {"TYR"," OH ",hyd_O,POL},
    {"TYR"," X  ",r_unk,UNK},
// 180
    {"VAL"," C  ",car_C,POL}, {"VAL"," O  ",car_O,POL}, {"VAL"," CA ",ali_C,APO},
    {"VAL"," N  ",ami_N,POL}, {"VAL"," CB ",ali_C,APO}, {"VAL"," CG1",ali_C,APO},
    {"VAL"," CG2",ali_C,APO}, {"VAL"," X  ",r_unk,UNK},
// 188
    {"ASX"," C  ",car_C,POL}, {"ASX"," O  ",car_O,POL}, {"ASX"," CA ",ali_C,APO},
    {"ASX"," N  ",ami_N,POL}, {"ASX"," CB ",ali_C,APO}, {"ASX"," CG ",car_C,POL},
    {"ASX"," XD1",unk_P,POL}, {"ASX"," XD2",unk_P,POL}, {"ASX"," AD1",unk_P,POL},
    {"ASX"," AD2",unk_P,POL}, {"ASX"," Y  ",r_unk,UNK},
// 199
    {"GLX"," C  ",car_C,POL}, {"GLX"," O  ",car_O,POL}, {"GLX"," CA ",ali_C,APO},
    {"GLX"," N  ",ami_N,POL}, {"GLX"," CB ",ali_C,APO}, {"GLX"," CG ",ali_C,APO},
    {"GLX"," CD ",car_C,POL}, {"GLX"," XE1",unk_P,POL}, {"GLX"," XE2",unk_P,POL},
    {"GLX"," AE1",unk_P,POL}, {"GLX"," AE2",unk_P,POL}, {"GLX"," Y  ",r_unk,UNK},
// 211
    {"CSE"," C  ",car_C,POL}, {"CSE"," O  ",car_O,POL}, {"CSE"," CA ",ali_C,APO},
    {"CSE"," N  ",ami_N,POL}, {"CSE"," CB ",ali_C,APO}, {"CSE"," SE ",sele,POL},
    {"CSE"," X  ",r_unk,UNK},
// 218
};

// tests freesasa_classify_radius() and freesasa_classify_oons_radius()
START_TEST (test_radius)
{
    char buf[50];
    for (int i = 0; i < n_atom_types; ++i) {
        const struct atom a = atoms[i];
        int oons_type = freesasa_classify_oons(a.a,a.b);
        double r1 = freesasa_classify_radius(a.a,a.b),
            r2 = freesasa_classify_oons_radius(oons_type);
        sprintf(buf,"%s %s ret=%f ref=%f",a.a,a.b,r1,a.radius);
        // make sure correct radius is supplied
        ck_assert_msg(fabs(r1 - a.radius) < 1e-10,buf);
        // make sure all regular atoms are given OONS-radii
        ck_assert(fabs(r1-r2) < 1e-10);
    }
    freesasa_set_verbosity(FREESASA_V_SILENT);
    // non OONS-atoms
    ck_assert(fabs(freesasa_classify_radius("XXX"," C  ") - 1.7) < 1e-10);
    ck_assert(fabs(freesasa_classify_radius("XXX"," N  ") - 1.55) < 1e-10 );
    ck_assert(fabs(freesasa_classify_radius("XXX"," S  ") - 1.8)  < 1e-10);
    ck_assert(fabs(freesasa_classify_radius("XXX"," O  ") - 1.52)  < 1e-10);
    ck_assert(fabs(freesasa_classify_radius("XXX"," P  ") - 1.8)  < 1e-10);
    ck_assert(fabs(freesasa_classify_radius("XXX"," X  ") - 0.0)  < 1e-10);
    //ck_assert(fabs(freesasa_classify_radius("XXX"," SE ") - 1.9)  < 1e-10);
    //ck_assert(fabs(freesasa_classify_radius("XXX"," H  ") - 1.2)  < 1e-10);
    freesasa_set_verbosity(0);
}
END_TEST

// tests freesasa_classify_class(), freesasa_classify_class2str(), ...oons2class(),
// ...noons() and ...nclasses(), ...oons2string()
START_TEST (test_class)
{
    char buf[50];
    for (int i = 0; i < n_atom_types; ++i) {
        const struct atom a = atoms[i];
        int c = freesasa_classify_class(a.a,a.b),
            o = freesasa_classify_oons(a.a,a.b);
        sprintf(buf,"%s %s %s",a.a,a.b,freesasa_classify_class2str(c));
        ck_assert_msg(freesasa_classify_class(a.a,a.b) == a.class, buf);
        ck_assert_msg(freesasa_classify_oons2class(o) == a.class,buf);

    }
    ck_assert(freesasa_classify_class("  C"," C1 ") == FREESASA_NUCLEICACID);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert(freesasa_classify_class("ABC"," X  ") == FREESASA_CLASS_UNKNOWN);
    freesasa_set_verbosity(0);
    ck_assert_str_eq(freesasa_classify_class2str(FREESASA_POLAR),"Polar");
    ck_assert_str_eq(freesasa_classify_class2str(FREESASA_APOLAR),"Apolar");
    ck_assert_str_eq(freesasa_classify_class2str(FREESASA_NUCLEICACID),"Nucleic");
    ck_assert_str_eq(freesasa_classify_class2str(FREESASA_CLASS_UNKNOWN),"Unknown");

    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert(freesasa_classify_class2str(100) == NULL);
    ck_assert(freesasa_classify_class2str(-1) == NULL);
    freesasa_set_verbosity(0);

    for (int i = 0; i < freesasa_classify_noons(); ++i) {
        ck_assert(freesasa_classify_oons2class(i) < freesasa_classify_nclasses());
    }
}
END_TEST

START_TEST (test_oons2str)
{

    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("ALA"," O  ")),"carbo_O");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("ALA"," C  ")),"carbo_C");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("ALA"," N  ")),"amide_N");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("ALA"," CA ")),"aliphatic_C");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("PHE"," CG ")),"aromatic_C");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("THR"," OG1")),"hydroxyl_O");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("CYS"," SG ")),"sulfur");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("CSE"," SE ")),"selenium");
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("ASX"," XD ")),"unknown_polar");
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_str_eq(freesasa_classify_oons2str(
                         freesasa_classify_oons("ABC"," X  ")),"unknown");
    ck_assert(freesasa_classify_oons2str(freesasa_classify_noons())==NULL);
    ck_assert(freesasa_classify_oons2str(-1)==NULL);
    freesasa_set_verbosity(0);

}
END_TEST

// tests freesasa_classify_residue2str(), freesasa_classify_residue(),
// freesasa_is_aminoacid() and freesasa_is_nucleicacid()
START_TEST (test_residue)
{
    // check consistency of classification and string output
    int nrt = freesasa_classify_nresiduetypes();
    const char **res = (const char**) malloc(nrt*sizeof(char*));
    for (int i = 0; i < nrt; ++i) {
        res[i] = freesasa_classify_residue2str(i);
        ck_assert_int_eq(freesasa_classify_residue(res[i]),i);
    }

    freesasa_set_verbosity(FREESASA_V_SILENT);
    // check erroneous input
    ck_assert(freesasa_classify_residue2str(nrt+1) == NULL);
    ck_assert(freesasa_classify_residue2str(-1) == NULL);
    ck_assert(freesasa_classify_is_aminoacid(-1) == FREESASA_FAIL);
    ck_assert(freesasa_classify_is_nucleicacid(-1) == FREESASA_FAIL);
    ck_assert(freesasa_classify_is_aminoacid(freesasa_classify_nresiduetypes())
              == FREESASA_FAIL);
    ck_assert(freesasa_classify_is_nucleicacid(freesasa_classify_nresiduetypes())
              == FREESASA_FAIL);
    ck_assert(freesasa_classify_residue("AAAA") == freesasa_classify_residue("UNK"));
    freesasa_set_verbosity(0);

    // check numbering
    int c;
    for (int i = 0; i < naa; ++i) {
        ck_assert((c = freesasa_classify_residue(aa[i])) < naa);
        ck_assert(freesasa_classify_is_aminoacid(c));
        ck_assert(!freesasa_classify_is_nucleicacid(c));
    }

    // irregular entries
    ck_assert((c = freesasa_classify_residue("UNK")) >= naa);
    ck_assert(!freesasa_classify_is_aminoacid(c));
    ck_assert(!freesasa_classify_is_nucleicacid(c));
    for (int i = 0; i < n_other; ++i) {
        ck_assert((c = freesasa_classify_residue(other_aa[i])) >= naa);
        ck_assert(freesasa_classify_is_aminoacid(c));
        ck_assert(!freesasa_classify_is_nucleicacid(c));
    }

    // nucleic acids
    for (int i = 0; i < n_nuc; ++i) {
        ck_assert((c = freesasa_classify_residue(nuc[i])) >= naa);
        ck_assert(!freesasa_classify_is_aminoacid(c));
        ck_assert(freesasa_classify_is_nucleicacid(c));
    }
}
END_TEST

// tests freesasa_classify_element() and freesasa_classify_element2str()
START_TEST (test_element)
{
    int c = freesasa_classify_element(" C ");
    ck_assert(c >= 0 && c < freesasa_classify_nelements());
    c = freesasa_classify_element(" N  ");
    ck_assert(c >= 0 && c < freesasa_classify_nelements());
    c = freesasa_classify_element(" O  ");
    ck_assert(c >= 0 && c < freesasa_classify_nelements());
    c = freesasa_classify_element(" S  ");
    ck_assert(c >= 0 && c < freesasa_classify_nelements());
    c = freesasa_classify_element(" P  ");
    ck_assert(c >= 0 && c < freesasa_classify_nelements());

    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" C  ")),"C");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" CD ")),"C");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" CE2")),"C");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element("   C    ")),"C");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element("   CE3  ")),"C");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" N  ")),"N");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" ND ")),"N");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" NE2")),"N");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" O  ")),"O");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" OD ")),"O");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" OE2")),"O");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" S  ")),"S");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" SD ")),"S");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" SE2")),"S");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" P  ")),"P");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" PD ")),"P");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element(" PE2")),"P");
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element("XXXX")),"unknown");
    ck_assert_str_eq(freesasa_classify_element2str(freesasa_classify_element("XXXXX")),"unknown");
    ck_assert(freesasa_classify_element2str(-1)==NULL);
    freesasa_set_verbosity(0);
}
END_TEST

START_TEST (test_user)
{
    FILE *f = fopen(DATADIR "oons.config","r");
    freesasa_classify* c = freesasa_classify_user(f);
    fclose(f);
    ck_assert(c != NULL);
    ck_assert(freesasa_classify_user_n_classes(c) == 2);
    ck_assert(fabs(freesasa_classify_user_radius(c,"ALA","CA") - 2.0) < 1e-5);
    ck_assert(fabs(freesasa_classify_user_radius(c,"ALA","N") - 1.55) < 1e-5);
    ck_assert_str_eq(freesasa_classify_user_class2str(c,freesasa_classify_user_class(c,"ALA","CB")), "apolar");
    ck_assert_str_eq(freesasa_classify_user_class2str(c,freesasa_classify_user_class(c,"ALA","O")), "polar");
    // compare oons.config and built in classification (should be identical for standard atoms)
    for (int i = 0; i < 188; ++i) {
        const char *res_name = atoms[i].a, *atom_name = atoms[i].b;
        if (strcmp(atom_name," X  ") == 0) continue;
        if (strcmp(atom_name," Y  ") == 0) continue;
        ck_assert(fabs(freesasa_classify_user_radius(c,res_name,atom_name) -
                       freesasa_classify_radius(res_name,atom_name)) < 1e-5);
        char *c1 = strdup(freesasa_classify_class2str(freesasa_classify_class(res_name,atom_name)));
        char *c2 = strdup(freesasa_classify_user_class2str(c,freesasa_classify_user_class(c,res_name,atom_name)));
        for (int i = 0; c1[i]; ++i) c1[i] = tolower(c1[i]); 
        for (int i = 0; c2[i]; ++i) c2[i] = tolower(c2[i]); 
        ck_assert_str_eq(c1,c2);
    }
    freesasa_classify_user_free(c);
}
END_TEST

Suite* classify_suite()
{
    Suite *s = suite_create("Classify");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core,test_radius);
    tcase_add_test(tc_core,test_class);
    tcase_add_test(tc_core,test_oons2str);
    tcase_add_test(tc_core,test_residue);
    tcase_add_test(tc_core,test_element);
    tcase_add_test(tc_core,test_user);
    suite_add_tcase(s,tc_core);

    return s;
}
