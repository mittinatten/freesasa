#include <stdio.h>
#include <stdlib.h>
#include "src/protein.h"
#include "src/sasa.h"

/* For now this program just outputs some raw data for sanity
   check. */

int main (int argc, char **argv) {
    protein *p = protein_init_from_pdb(stdin);
    double *sasa = (double*) malloc(sizeof(double)*protein_n(p));
    double *r = (double*) malloc(sizeof(double)*protein_n(p));
    protein_r_def(p,r);
    sasa_shrake_rupley(sasa,protein_xyz(p),r,protein_n(p),1000);
    printf ("%ul\n",protein_n(p));
    for (int i = 0; i < protein_n(p); ++i) {
	printf("%5i,%4s,%3s,%4s,%c,%6.2f\n",
	       i,
	       protein_atom_name(p,i),
	       protein_atom_res_name(p,i),
	       protein_atom_res_number(p,i),
	       protein_atom_chain(p,i),
	       //atom_type2str(a->at),
	       //atom_class2str(a->ac),
	       sasa[i]);
	
    }
    protein_free(p);
    free(sasa);
    free(r);
}
