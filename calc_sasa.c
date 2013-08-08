/*
  Copyright Simon Mitternacht 2013.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/protein.h"
#include "src/sasa.h"
#include "src/oons.h"
#include "src/srp.h"

#define DEF_SR_POINTS 100
#define DEF_LR_SPACING 0.25

enum algorithms {LEE_RICHARDS, SHRAKE_RUPLEY};
const char *alg_names[] = {"Lee & Richards", "Shrake & Rupley"};

#if __STDC__
extern int getopt(int, char * const *, const char *);
extern int optind;
#endif

void help(const char* argv0) {
    fprintf(stderr,"\nUsage: %s [options] pdb-file(s)\n",argv0);
    fprintf(stderr,"\nOptions are:\n"
	    "       -h  print this message\n"
	    "       -S  use Shrake & Rupley alogrithm [default]\n"
	    "       -n  number of test points in Shrake & Rupley algorithm\n"
	    "           Default is %d, allowed values are:\n"
	    "           ",DEF_SR_POINTS);
    srp_print_n_opt(stderr);
    fprintf(stderr,"       -L  use Lee & Richards algorithm\n"
	           "       -d  grid spacing in Lee & Richards algorithm\n"
	           "           Default value is %4.2f Å\n"
	    ,DEF_LR_SPACING);
    fprintf(stderr,"\nIf no pdb-file is specified STDIN is used for input.\n\n");
}

void short_help(const char* argv0) {
    fprintf(stderr,"Run '%s -h' for help.\n\n", argv0);
}

void run_analysis(FILE *input, int use_alg, const char *name, void *param) {
    protein *p; 
    double *sasa, *r;
    p  = protein_init_from_pdb(input);
    r = (double*) malloc(sizeof(double)*protein_n(p));
    sasa = (double*) malloc(sizeof(double)*protein_n(p));

    printf("# Using van der Waals radii and atom classes defined \n"
	   "# by Ooi et al (PNAS 1987, 84:3086-3090) and a probe raidus\n"
	   "# of %f Å.\n\n", PROBE_RADIUS);

    printf("File: %s\n",name);

    //calc OONS radii
    protein_r_def(p,r);

    printf("algorithm: %s\n",alg_names[use_alg]);

    switch(use_alg) {
    case SHRAKE_RUPLEY:
	printf("N_testpoint: %d\n",*(int *)param);
	sasa_shrake_rupley(sasa,protein_xyz(p),r,protein_n(p),*(int *)param);
	break;
    case LEE_RICHARDS:
	printf("d_slice: %f Å.\n",*(double *)param);
	sasa_lee_richards(sasa,protein_xyz(p),r,protein_n(p),*(double *)param);
	break;
    default:
	fprintf(stderr,"Error: no SASA algorithm specified.\n");
	exit(0);
    }
    sasa_per_atomclass(stdout,oons_classes(),p,sasa);
    sasa_per_atomclass(stdout,oons_types(),p,sasa);
    sasa_per_atomclass(stdout,atomclassifier_residue(),p,sasa);
    
    protein_free(p);
    free(sasa);
    free(r);    
}

int main (int argc, char **argv) {

    int n_sr_points = DEF_SR_POINTS;
    double d_lr = DEF_LR_SPACING;
    void *param;
    int use_alg = SHRAKE_RUPLEY, alg_set = 0;
    FILE *input = NULL;

    extern char *optarg;
    char opt;

    while ((opt = getopt(argc, argv, "f:n:d:hLS")) != -1) {
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
	    break;
	case 'S':
	    use_alg = SHRAKE_RUPLEY;
	    ++alg_set;
	    break;
	case 'L':
	    use_alg = LEE_RICHARDS;
	    ++alg_set;
	    break;
	case 'd':
	    d_lr = atof(optarg);
	    break;
	default:
	    fprintf(stderr, "\nWarning: unknown option '%c' (will be ignored)\n\n", 
		    opt);
	    break;
	}
    }
    if (alg_set > 1) {
	fprintf(stderr, "Error: multiple algorithms specified.\n");
	exit(0);
    }
    if (use_alg == LEE_RICHARDS) param = &d_lr;
    if (use_alg == SHRAKE_RUPLEY) param = &n_sr_points;

    printf("\n# SASALIB 2013\n");
    
    if (argc > 0) {
	for (int i = optind; i < argc; ++i) {
	    input = fopen(argv[i],"r");
	    if (input != NULL) {
		run_analysis(input,use_alg,argv[i],param);
		fclose(input);
	    } else {
		fprintf(stderr, "Error opening file '%s'\n", argv[i]);
		exit(0);
	    }	    
	}
    } else {
	run_analysis(stdin,use_alg,"stdin",param);
    }    
    return 0;
}
