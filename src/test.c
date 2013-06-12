#include <stdio.h>
#include <stdlib.h>
#include "protein.h"
#include "sasa.h"

/* For now this program just outputs some raw data for sanity
   check. */

int main (int argc, char **argv) {
    protein *p = protein_init_from_pdb(stdin);
    double *sasa = (double*) malloc(sizeof(double)*p->number_atoms);
    sasa_shrake_rupley(sasa,p->xyz,p->r,p->number_atoms,2000);
    for (int i = 0; i < p->number_atoms; ++i) {
	atom *a = &p->a[i];
	printf("%5i,%4s,%3s,%4s,%c,%12s,%10s,%6.2f\n",
	       i,a->atom_name,a->res_name,
	       a->res_number,a->chain_label,
	       atom_type2str(a->at),
	       atom_class2str(a->ac),
	       sasa[i]);
	
    }
    protein_free(p);
    free(sasa);
}
