#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/protein.h"
#include "src/sasa.h"

#define DEF_SR_POINTS 100

#if __STDC__
extern int getopt(int, char * const *, const char *);
#endif

int main (int argc, char **argv) {
    int n_sr_points = DEF_SR_POINTS;
    protein *p; 
    double *sasa;
    FILE *input = stdin;

    extern char *optarg;
    char opt;

    while ((opt = getopt(argc, argv, "f:n:")) != -1) {
        switch(opt) {
	case 'f':
	    input = fopen(optarg, "r");
	    if (input == NULL) {
		fprintf(stderr,"\nError: could not open file '%s'.\n", 
			optarg);
		exit(1);
	    }
	    break;
	case 'n':
	    n_sr_points = atoi(optarg);
	    if (n_sr_points > MAX_SR_POINTS) {
		fprintf(stderr,"\nError: n = %d is too large (%d is limit).\n",
			n_sr_points, MAX_SR_POINTS);
		exit(1);
	    }
	    break;
	default:
	    fprintf(stderr, "\nWarning: unknown option '%c' (will be ignored)\n\n", 
		    opt);
	    break;
	}
    }

    p  = protein_init_from_pdb(input);
    sasa = (double*) malloc(sizeof(double)*p->number_atoms);

    sasa_shrake_rupley(sasa,p->xyz,p->r,p->number_atoms,n_sr_points);

    protein_free(p);
    free(sasa);
    fclose(input);
}
