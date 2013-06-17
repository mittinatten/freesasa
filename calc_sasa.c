#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/protein.h"
#include "src/sasa.h"

#define DEF_SR_POINTS 100

#if __STDC__
extern int getopt(int, char * const *, const char *);
#endif

void help(const char* argv0) {
    fprintf(stderr,"\nUsage: %s [options] [< pdb-file]\n",argv0);
    fprintf(stderr,"\nOptions are:\n"
	           "       -h  print this message\n"
                   "       -f  pdb-file\n"
	           "       -n  number of test points in Shrake & Rupley algorithm\n"
	           "           Default is %d, max value is %d.\n"
	    ,DEF_SR_POINTS,MAX_SR_POINTS);
    fprintf(stderr,"\nIf no pdb-file is specified STDIN is used for input.\n\n");
}

void short_help(const char* argv0) {
    fprintf(stderr,"Run '%s -h' for help.\n\n", argv0);
}

int main (int argc, char **argv) {
    int n_sr_points = DEF_SR_POINTS;
    protein *p; 
    double *sasa;
    FILE *input = stdin;

    extern char *optarg;
    char opt;

    while ((opt = getopt(argc, argv, "f:n:h")) != -1) {
        switch(opt) {
	case 'h':
	    help(argv[0]);
	    exit(0);
	case 'f':
	    input = fopen(optarg, "r");
	    if (input == NULL) {
		fprintf(stderr,"\nError: could not open file '%s'.\n\n", 
			optarg);
		short_help(argv[0]);
		exit(1);
	    }
	    break;
	case 'n':
	    n_sr_points = atoi(optarg);
	    if (n_sr_points > MAX_SR_POINTS) {
		fprintf(stderr,"\nError: n = %d is too large (%d is limit).\n\n",
			n_sr_points, MAX_SR_POINTS);
		short_help(argv[0]);
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
    double *r = (double*) malloc(sizeof(double)*protein_n(p));
    sasa = (double*) malloc(sizeof(double)*protein_n(p));
    
    //calc OONS radii
    protein_r_def(p,r);

    sasa_shrake_rupley(sasa,protein_xyz(p),
		       r,protein_n(p),n_sr_points);

    protein_free(p);
    free(sasa);
    free(r);
    fclose(input);
}
