#include <string.h>
#include "integrate.h"

void integrate_sasa_per_atomclass(FILE *out, atomclassifier ac,
                                  protein *p, double *sasa)
{
    int nc = ac.nclasses;
    double s[nc];
    for (int i = 0; i < nc; ++i) {
        s[i] = 0;
    }
    for (size_t i = 0; i < protein_n(p); ++i) {
        int class = ac.classify(protein_atom_res_name(p,i),
                                protein_atom_name(p,i));
        s[class] += sasa[i];
    }
    for (int i = 0; i < nc; ++i) {
        fprintf(out,"%s\t%6.2f\n",ac.class2str[i],s[i]);
    }
}


