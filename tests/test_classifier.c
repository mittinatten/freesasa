#include "tools.h"
#include <check.h>
#include <classifier.h>
#include <freesasa.h>
#include <freesasa_internal.h>
#include <stdio.h>
#include <stdlib.h>

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
#define phos 1.80
#define unk_P 1.5
#define r_unk -1.0
#define APO FREESASA_ATOM_APOLAR
#define POL FREESASA_ATOM_POLAR
#define UNK FREESASA_ATOM_UNKNOWN
#define NUC 3

#define naa 20
const char *aa[naa] = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                       "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                       "PRO", "SER", "THR", "TRP", "TYR", "VAL"};
#define n_other 7
const char *other_aa[n_other] = {"GLX", "SEC", "MSE", "PYL", "PYH", "ASX", "CSE"};

#define n_capping 2
const char *capping[n_capping] = {"ACE", "NH2"};

#define n_nuc 13
const char *nuc[n_nuc] = {"DA", "DC", "DG", "DT", "DU", "DI", "A", "C", "G", "T", "U", "I", "N"};

// clang-format off
const struct atom atoms[] = {
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
    {"MET"," SD ",sulf, POL}, {"MET"," CE ",ali_C,APO}, {"MET"," X  ",r_unk,UNK},
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
    {"CSE"," N  ",ami_N,POL}, {"CSE"," CB ",ali_C,APO}, {"CSE","SE  ",sele,POL},
    {"CSE"," X  ",r_unk,UNK},
// 218
    {"SEC"," C  ",car_C,POL}, {"SEC"," O  ",car_O,POL}, {"SEC"," CA ",ali_C,APO},
    {"SEC"," N  ",ami_N,POL}, {"SEC"," CB ",ali_C,APO}, {"SEC","SE  ",sele,POL},
    {"SEC"," X  ",r_unk,UNK},
// 225
    {"PYL"," C  ",car_C,POL}, {"PYL"," O  ",car_O,POL}, {"PYL"," CA ",ali_C,APO},
    {"PYL"," N  ",ami_N,POL}, {"PYL"," CB ",ali_C,APO}, {"PYL"," CG ",ali_C,APO},
    {"PYL"," CG ",ali_C,APO}, {"PYL"," CD ",ali_C,APO}, {"PYL"," CE ",ali_C,APO},
    {"PYL"," NZ ",ami_N,POL}, {"PYL"," O2 ",car_O,POL}, {"PYL"," C2 ",car_C,POL},
    {"PYL"," CA2",aro_C,APO}, {"PYL"," CB2",ali_C,APO}, {"PYL"," CG2",aro_C,APO},
    {"PYL"," CD2",aro_C,APO}, {"PYL"," CE2",aro_C,APO}, {"PYL"," N2 ",ami_N,POL},
// 243
    {"ACE"," C  ",car_C,POL}, {"ACE"," O  ",car_O,POL}, {"ACE"," CH3",ali_C,APO},
    {"NH2"," NH2",ami_N,POL},
// 247
    {"  A"," P  ",phos, POL}, {"  A"," OP1",car_O,POL}, {"  A"," OP2",car_O,POL},
    {"  A"," OP3",car_O,POL}, {"  A"," O5'",car_O,POL}, {"  A"," C5'",ali_C,APO},
    {"  A"," C4'",aro_C,APO}, {"  A"," O4'",car_O,POL}, {"  A"," C3'",aro_C,APO},
    {"  A"," O3'",car_O,POL}, {"  A"," C2'",aro_C,APO}, {"  A"," O2'",car_O,POL},
    {"  A"," C1'",aro_C,APO}, {"  A"," N9 ",ami_N,POL}, {"  A"," C8 ",aro_C,APO},
    {"  A"," N7 ",ami_N,POL}, {"  A"," C5 ",aro_C,APO}, {"  A"," C6 ",aro_C,APO},
    {"  A"," N6 ",ami_N,POL}, {"  A"," N1 ",ami_N,POL}, {"  A"," C2 ",aro_C,APO},
    {"  A"," N3 ",ami_N,POL}, {"  A"," C4 ",aro_C,APO},
// 270
    {"  C"," P  ",phos, POL}, {"  C"," OP1",car_O,POL}, {"  C"," OP2",car_O,POL},
    {"  C"," OP3",car_O,POL}, {"  C"," O5'",car_O,POL}, {"  C"," C5'",ali_C,APO},
    {"  C"," C4'",aro_C,APO}, {"  C"," O4'",car_O,POL}, {"  C"," C3'",aro_C,APO},
    {"  C"," O3'",car_O,POL}, {"  C"," C2'",aro_C,APO}, {"  C"," O2'",car_O,POL},
    {"  C"," C1'",aro_C,APO}, {"  C"," N1 ",ami_N,POL}, {"  C"," C2 ",aro_C,APO},
    {"  C"," O2 ",car_O,POL}, {"  C"," N3 ",ami_N,POL}, {"  C"," C4 ",aro_C,APO},
    {"  C"," N4 ",ami_N,POL}, {"  C"," C5 ",aro_C,APO}, {"  C"," C6 ",aro_C,APO},
// 291
    {"  G"," P  ",phos, POL}, {"  G"," OP1",car_O,POL}, {"  G"," OP2",car_O,POL},
    {"  G"," OP3",car_O,POL}, {"  G"," O5'",car_O,POL}, {"  G"," C5'",ali_C,APO},
    {"  G"," C4'",aro_C,APO}, {"  G"," O4'",car_O,POL}, {"  G"," C3'",aro_C,APO},
    {"  G"," O3'",car_O,POL}, {"  G"," C2'",aro_C,APO}, {"  G"," O2'",car_O,POL},
    {"  G"," C1'",aro_C,APO}, {"  G"," N9 ",ami_N,POL}, {"  G"," C8 ",aro_C,APO},
    {"  G"," N7 ",ami_N,POL}, {"  G"," C5 ",aro_C,APO}, {"  G"," C6 ",aro_C,APO},
    {"  G"," O6 ",car_O,POL}, {"  G"," N1 ",ami_N,POL}, {"  G"," C2 ",aro_C,APO},
    {"  G"," N3 ",ami_N,POL}, {"  G"," C4 ",aro_C,APO},
// 304
    {"  T"," P  ",phos, POL}, {"  T"," OP1",car_O,POL}, {"  T"," OP2",car_O,POL},
    {"  T"," OP3",car_O,POL}, {"  T"," O5'",car_O,POL}, {"  T"," C5'",ali_C,APO},
    {"  T"," C4'",aro_C,APO}, {"  T"," O4'",car_O,POL}, {"  T"," C3'",aro_C,APO},
    {"  T"," O3'",car_O,POL}, {"  T"," C2'",aro_C,APO}, {"  T"," O2'",car_O,POL},
    {"  T"," C1'",aro_C,APO}, {"  T"," N1 ",ami_N,POL}, {"  T"," C2 ",aro_C,APO},
    {"  T"," O2 ",car_O,POL}, {"  T"," N3 ",ami_N,POL}, {"  T"," C4 ",aro_C,APO},
    {"  T"," O4 ",car_O,POL}, {"  T"," C5 ",aro_C,APO}, {"  T"," C6 ",aro_C,APO},
    {"  T"," C7 ",aro_C,APO},
// 326
    {"  U"," P  ",phos, POL}, {"  U"," OP1",car_O,POL}, {"  U"," OP2",car_O,POL},
    {"  U"," OP3",car_O,POL}, {"  U"," O5'",car_O,POL}, {"  U"," C5'",ali_C,APO},
    {"  U"," C4'",aro_C,APO}, {"  U"," O4'",car_O,POL}, {"  U"," C3'",aro_C,APO},
    {"  U"," O3'",car_O,POL}, {"  U"," C2'",aro_C,APO}, {"  U"," O2'",car_O,POL},
    {"  U"," C1'",aro_C,APO}, {"  U"," N1 ",ami_N,POL}, {"  U"," C2 ",aro_C,APO},
    {"  U"," O2 ",car_O,POL}, {"  U"," N3 ",ami_N,POL}, {"  U"," C4 ",aro_C,APO},
    {"  U"," O4 ",car_O,POL}, {"  U"," C5 ",aro_C,APO}, {"  U"," C6 ",aro_C,APO},
// 347
    {"  I"," P  ",phos, POL}, {"  I"," OP1",car_O,POL}, {"  I"," OP2",car_O,POL},
    {"  I"," OP3",car_O,POL}, {"  I"," O5'",car_O,POL}, {"  I"," C5'",ali_C,APO},
    {"  I"," C4'",aro_C,APO}, {"  I"," O4'",car_O,POL}, {"  I"," C3'",aro_C,APO},
    {"  I"," O3'",car_O,POL}, {"  I"," C2'",aro_C,APO}, {"  I"," O2'",car_O,POL},
    {"  I"," C1'",aro_C,APO}, {"  I"," N9 ",ami_N,POL}, {"  I"," C8 ",aro_C,APO},
    {"  I"," N7 ",ami_N,POL}, {"  I"," C5 ",aro_C,APO}, {"  I"," C6 ",aro_C,APO},
    {"  I"," N6 ",ami_N,POL}, {"  I"," N1 ",ami_N,POL}, {"  I"," C2 ",aro_C,APO},
    {"  I"," N3 ",ami_N,POL}, {"  I"," C4 ",aro_C,APO},

    {" DA"," P  ",phos, POL}, {" DA"," OP1",car_O,POL}, {" DA"," OP2",car_O,POL},
    {" DA"," OP3",car_O,POL}, {" DA"," O5'",car_O,POL}, {" DA"," C5'",ali_C,APO},
    {" DA"," C4'",aro_C,APO}, {" DA"," O4'",car_O,POL}, {" DA"," C3'",aro_C,APO},
    {" DA"," O3'",car_O,POL}, {" DA"," C2'",aro_C,APO},
    {" DA"," C1'",aro_C,APO}, {" DA"," N9 ",ami_N,POL}, {" DA"," C8 ",aro_C,APO},
    {" DA"," N7 ",ami_N,POL}, {" DA"," C5 ",aro_C,APO}, {" DA"," C6 ",aro_C,APO},
    {" DA"," N6 ",ami_N,POL}, {" DA"," N1 ",ami_N,POL}, {" DA"," C2 ",aro_C,APO},
    {" DA"," N3 ",ami_N,POL}, {" DA"," C4 ",aro_C,APO},

    {" DC"," P  ",phos, POL}, {" DC"," OP1",car_O,POL}, {" DC"," OP2",car_O,POL},
    {" DC"," OP3",car_O,POL}, {" DC"," O5'",car_O,POL}, {" DC"," C5'",ali_C,APO},
    {" DC"," C4'",aro_C,APO}, {" DC"," O4'",car_O,POL}, {" DC"," C3'",aro_C,APO},
    {" DC"," O3'",car_O,POL}, {" DC"," C2'",aro_C,APO},
    {" DC"," C1'",aro_C,APO}, {" DC"," N1 ",ami_N,POL}, {" DC"," C2 ",aro_C,APO},
    {" DC"," O2 ",car_O,POL}, {" DC"," N3 ",ami_N,POL}, {" DC"," C4 ",aro_C,APO},
    {" DC"," N4 ",ami_N,POL}, {" DC"," C5 ",aro_C,APO}, {" DC"," C6 ",aro_C,APO},

    {" DG"," P  ",phos, POL}, {" DG"," OP1",car_O,POL}, {" DG"," OP2",car_O,POL},
    {" DG"," OP3",car_O,POL}, {" DG"," O5'",car_O,POL}, {" DG"," C5'",ali_C,APO},
    {" DG"," C4'",aro_C,APO}, {" DG"," O4'",car_O,POL}, {" DG"," C3'",aro_C,APO},
    {" DG"," O3'",car_O,POL}, {" DG"," C2'",aro_C,APO},
    {" DG"," C1'",aro_C,APO}, {" DG"," N9 ",ami_N,POL}, {" DG"," C8 ",aro_C,APO},
    {" DG"," N7 ",ami_N,POL}, {" DG"," C5 ",aro_C,APO}, {" DG"," C6 ",aro_C,APO},
    {" DG"," O6 ",car_O,POL}, {" DG"," N1 ",ami_N,POL}, {" DG"," C2 ",aro_C,APO},
    {" DG"," N3 ",ami_N,POL}, {" DG"," C4 ",aro_C,APO},

    {" DT"," P  ",phos, POL}, {" DT"," OP1",car_O,POL}, {" DT"," OP2",car_O,POL},
    {" DT"," OP3",car_O,POL}, {" DT"," O5'",car_O,POL}, {" DT"," C5'",ali_C,APO},
    {" DT"," C4'",aro_C,APO}, {" DT"," O4'",car_O,POL}, {" DT"," C3'",aro_C,APO},
    {" DT"," O3'",car_O,POL}, {" DT"," C2'",aro_C,APO},
    {" DT"," C1'",aro_C,APO}, {" DT"," N1 ",ami_N,POL}, {" DT"," C2 ",aro_C,APO},
    {" DT"," O2 ",car_O,POL}, {" DT"," N3 ",ami_N,POL}, {" DT"," C4 ",aro_C,APO},
    {" DT"," O4 ",car_O,POL}, {" DT"," C5 ",aro_C,APO}, {" DT"," C6 ",aro_C,APO},
    {" DT"," C7 ",aro_C,APO},

    {" DU"," P  ",phos, POL}, {" DU"," OP1",car_O,POL}, {" DU"," OP2",car_O,POL},
    {" DU"," OP3",car_O,POL}, {" DU"," O5'",car_O,POL}, {" DU"," C5'",ali_C,APO},
    {" DU"," C4'",aro_C,APO}, {" DU"," O4'",car_O,POL}, {" DU"," C3'",aro_C,APO},
    {" DU"," O3'",car_O,POL}, {" DU"," C2'",aro_C,APO},
    {" DU"," C1'",aro_C,APO}, {" DU"," N1 ",ami_N,POL}, {" DU"," C2 ",aro_C,APO},
    {" DU"," O2 ",car_O,POL}, {" DU"," N3 ",ami_N,POL}, {" DU"," C4 ",aro_C,APO},
    {" DU"," O4 ",car_O,POL}, {" DU"," C5 ",aro_C,APO}, {" DU"," C6 ",aro_C,APO},

    {" DI"," P  ",phos, POL}, {" DI"," OP1",car_O,POL}, {" DI"," OP2",car_O,POL},
    {" DI"," OP3",car_O,POL}, {" DI"," O5'",car_O,POL}, {" DI"," C5'",ali_C,APO},
    {" DI"," C4'",aro_C,APO}, {" DI"," O4'",car_O,POL}, {" DI"," C3'",aro_C,APO},
    {" DI"," O3'",car_O,POL}, {" DI"," C2'",aro_C,APO},
    {" DI"," C1'",aro_C,APO}, {" DI"," N9 ",ami_N,POL}, {" DI"," C8 ",aro_C,APO},
    {" DI"," N7 ",ami_N,POL}, {" DI"," C5 ",aro_C,APO}, {" DI"," C6 ",aro_C,APO},
    {" DI"," N6 ",ami_N,POL}, {" DI"," N1 ",ami_N,POL}, {" DI"," C2 ",aro_C,APO},
    {" DI"," N3 ",ami_N,POL}, {" DI"," C4 ",aro_C,APO},

    {"MSE"," C  ",car_C,POL}, {"MSE"," O  ",car_O,POL}, {"MSE"," CA ",ali_C,APO},
    {"MSE"," N  ",ami_N,POL}, {"MSE"," CB ",ali_C,APO}, {"MSE"," CG ",ali_C,APO},
    {"MSE"," SE ",sele, POL}, {"MSE"," CE ",ali_C,APO}, {"MSE"," X  ",r_unk,UNK},
};
// clang-format on

const int n_atom_types = sizeof(atoms) / sizeof(struct atom);

static const freesasa_classifier *oons_c = &freesasa_oons_classifier;

// tests freesasa_classifier_radius() and freesasa_guess_radius
START_TEST(test_radius)
{
    char buf[50];
    freesasa_set_verbosity(FREESASA_V_SILENT); // the X and Y atoms give warnings
    for (int i = 0; i < n_atom_types; ++i) {
        const struct atom a = atoms[i];
        double r1 = freesasa_classifier_radius(oons_c, a.a, a.b);
        ck_assert(float_eq(r1, a.radius, 1e-10));
    }
    ck_assert(float_eq(freesasa_guess_radius(" H"), 1.10, 1e-10));
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST(test_class)
{
    char buf[50];
    for (int i = 0; i < n_atom_types; ++i) {
        const struct atom a = atoms[i];
        int c = freesasa_classifier_class(oons_c, a.a, a.b);
        if (c == FREESASA_WARN) {
            ck_assert(a.class == UNK);
            continue;
        }
        sprintf(buf, "Classification error for %s %s %s %s",
                a.a, a.b, freesasa_classifier_class2str(c),
                freesasa_classifier_class2str(a.class));
        ck_assert_msg(c == a.class, buf);
    }
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert(freesasa_classifier_class(oons_c, "ABC", " X  ") == FREESASA_ATOM_UNKNOWN);

    ck_assert(freesasa_classifier_class2str(100) == NULL);
    ck_assert(freesasa_classifier_class2str(-1) == NULL);
    ck_assert_str_eq(freesasa_classifier_class2str(FREESASA_ATOM_APOLAR), "Apolar");
    ck_assert_str_eq(freesasa_classifier_class2str(FREESASA_ATOM_POLAR), "Polar");
    ck_assert_str_eq(freesasa_classifier_class2str(FREESASA_ATOM_UNKNOWN), "Unknown");
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST(test_residue)
{
    // check consistency of classification and string output
    int nrt = freesasa_classify_n_residue_types();
    for (int i = 0; i < nrt; ++i) {
        const char *res = freesasa_classify_residue_name(i);
        ck_assert_int_eq(freesasa_classify_residue(res), i);
    }

    // check numbering
    for (int i = 0; i < naa; ++i) {
        ck_assert_int_lt(freesasa_classify_residue(aa[i]), naa);
    }

    ck_assert_int_ge(freesasa_classify_residue("UNK"), naa);
    ck_assert_int_ge(freesasa_classify_residue("XXX"), naa);

    for (int i = 0; i < n_other; ++i) {
        ck_assert_int_ge(freesasa_classify_residue(other_aa[i]), naa);
    }
    for (int i = 0; i < n_capping; ++i) {
        ck_assert_int_ge(freesasa_classify_residue(capping[i]), naa);
    }
    for (int i = 0; i < n_nuc; ++i) {
        ck_assert_int_ge(freesasa_classify_residue(nuc[i]), naa);
    }
}
END_TEST

START_TEST(test_user)
{
    FILE *clf = fopen(SHAREDIR "oons.config", "r");
    ck_assert(clf != NULL);

    freesasa_classifier *c = freesasa_classifier_from_file(clf);
    ck_assert(c != NULL);
    ck_assert(freesasa_classifier_class(c, "ALA", "CA") == FREESASA_ATOM_APOLAR);
    ck_assert(freesasa_classifier_class(c, "ALA", "O") == FREESASA_ATOM_POLAR);
    ck_assert(float_eq(freesasa_classifier_radius(c, "ALA", "CA"), 2.0, 1e-5));
    ck_assert(float_eq(freesasa_classifier_radius(c, "ALA", "N"), 1.55, 1e-5));
    // compare oons.config and built in classification (should be identical for standard atoms)
    for (int i = 0; i < 188; ++i) {
        const char *res_name = atoms[i].a, *atom_name = atoms[i].b;
        if (strcmp(atom_name, " X  ") == 0) continue;
        if (strcmp(atom_name, " Y  ") == 0) continue;
        ck_assert(float_eq(freesasa_classifier_radius(c, res_name, atom_name),
                           freesasa_classifier_radius(oons_c, res_name, atom_name), 1e-5));
        ck_assert(freesasa_classifier_class(c, res_name, atom_name) ==
                  freesasa_classifier_class(oons_c, res_name, atom_name));
    }
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert(freesasa_classifier_radius(c, "ALA", "X") < 0);
    ck_assert(freesasa_classifier_radius(c, "X", "CB") > 0);
    ck_assert(freesasa_classifier_radius(c, "X", "X") < 0);
    ck_assert(freesasa_classifier_class(c, "ALA", "X") == FREESASA_ATOM_UNKNOWN);
    ck_assert(freesasa_classifier_class(c, "X", "CB") == FREESASA_ATOM_APOLAR);
    ck_assert(freesasa_classifier_class(c, "X", "X") == FREESASA_ATOM_UNKNOWN);
    freesasa_classifier_free(c);

    FILE *f = fopen(DATADIR "empty.pdb", "r");
    c = freesasa_classifier_from_file(f);
    ck_assert(c == NULL);
    fclose(f);
    freesasa_classifier_free(c);

    f = fopen(DATADIR "err.config", "r");
    c = freesasa_classifier_from_file(f);
    ck_assert(c == NULL);
    fclose(f);
    freesasa_classifier_free(c);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST(test_backbone)
{
    ck_assert(freesasa_atom_is_backbone("C"));
    ck_assert(freesasa_atom_is_backbone(" C  "));
    ck_assert(freesasa_atom_is_backbone(" CA "));
    ck_assert(freesasa_atom_is_backbone(" CA     "));
    ck_assert(freesasa_atom_is_backbone("CA"));
    ck_assert(freesasa_atom_is_backbone("N"));
    ck_assert(freesasa_atom_is_backbone("O"));
    ck_assert(freesasa_atom_is_backbone("OXT"));
    ck_assert(freesasa_atom_is_backbone("P"));
    ck_assert(freesasa_atom_is_backbone("OP1"));
    ck_assert(freesasa_atom_is_backbone("OP2"));
    ck_assert(freesasa_atom_is_backbone("O5'"));
    ck_assert(freesasa_atom_is_backbone("C5'"));
    ck_assert(freesasa_atom_is_backbone("C4'"));
    ck_assert(freesasa_atom_is_backbone("O4'"));
    ck_assert(freesasa_atom_is_backbone("C3'"));
    ck_assert(freesasa_atom_is_backbone("O3'"));
    ck_assert(freesasa_atom_is_backbone("C2'"));
    ck_assert(freesasa_atom_is_backbone("C1'"));
    ck_assert(freesasa_atom_is_backbone("") == 0);
    ck_assert(freesasa_atom_is_backbone("X") == 0);
}
END_TEST

START_TEST(test_memerr)
{
    freesasa_set_verbosity(FREESASA_V_SILENT);
    set_fail_after(1);
    void *ptr[] = {freesasa_classifier_types_new(),
                   freesasa_classifier_residue_new("ALA"),
                   freesasa_classifier_new()};
    set_fail_after(0);

    for (int i = 0; i < sizeof(ptr) / sizeof(void *); ++i) {
        ck_assert_ptr_eq(ptr[i], NULL);
    }

    for (int i = 1; i < 4; ++i) {
        struct classifier_types *types = freesasa_classifier_types_new();
        struct classifier_residue *res = freesasa_classifier_residue_new("ALA");
        struct freesasa_classifier *clf = freesasa_classifier_new();
        int ret;

        if (i < 2) {
            set_fail_after(i);
            ret = freesasa_classifier_add_type(types, "a", "A", 1.0);
            set_fail_after(0);
            ck_assert_int_eq(ret, FREESASA_FAIL);
        }
        set_fail_after(i);
        ret = freesasa_classifier_add_atom(res, "A", 1.0, 0);
        set_fail_after(0);
        ck_assert_int_eq(ret, FREESASA_FAIL);

        set_fail_after(i);
        ret = freesasa_classifier_add_residue(clf, "A");
        set_fail_after(0);
        ck_assert_int_eq(ret, FREESASA_FAIL);

        freesasa_classifier_types_free(types);
        freesasa_classifier_residue_free(res);
        freesasa_classifier_free(clf);
    }
    // don't test all levels, but make sure errors in low level
    // allocation propagates to the interface
    FILE *config = fopen(SHAREDIR "naccess.config", "r");
    for (int i = 1; i < 200; ++i) {
        set_fail_after(i);
        void *ptr = freesasa_classifier_from_file(config);
        set_fail_after(0);
        ck_assert_ptr_eq(ptr, NULL);
        rewind(config);
    }
    fclose(config);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

extern TCase *test_classifier_static();

Suite *classifier_suite()
{
    Suite *s = suite_create("Classify");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_radius);
    tcase_add_test(tc_core, test_class);
    tcase_add_test(tc_core, test_residue);
    tcase_add_test(tc_core, test_user);
    tcase_add_test(tc_core, test_backbone);
    tcase_add_test(tc_core, test_memerr);

    TCase *tc_static = test_classifier_static();

    suite_add_tcase(s, tc_core);
    suite_add_tcase(s, tc_static);

    return s;
}
