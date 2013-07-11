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
#endif

void help(const char* argv0) {
    fprintf(stderr,"\nUsage: %s [options] [< pdb-file]\n",argv0);
    fprintf(stderr,"\nOptions are:\n"
	    "       -h  print this message\n"
	    "       -f  pdb-file\n"
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

int main (int argc, char **argv) {

    int n_sr_points = DEF_SR_POINTS;
    double d_lr = DEF_LR_SPACING;
    int use_alg = SHRAKE_RUPLEY, alg_set = 0;
    protein *p; 
    double *sasa;
    FILE *input = stdin;

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
    p  = protein_init_from_pdb(input);
    double *r = (double*) malloc(sizeof(double)*protein_n(p));
    sasa = (double*) malloc(sizeof(double)*protein_n(p));
    
    printf("\n\n### SASALIB 2013 ###\n\n");
	   
    printf("# Using van der Waals radii and atom classes defined \n"
	   "# by Ooi et al (PNAS 1987, 84:3086-3090) and a probe raidus\n"
	   "# of %f Å.\n\n", PROBE_RADIUS);

    //calc OONS radii
    protein_r_def(p,r);

    switch(use_alg) {
    case SHRAKE_RUPLEY:
	printf("# Using Shrake & Rupley algorithm with %d test-points\n",n_sr_points);
	sasa_shrake_rupley(sasa,protein_xyz(p),r,protein_n(p),n_sr_points);
	break;
    case LEE_RICHARDS:
	printf("# Using Lee & Richards algorithm with grid spacing of %f Å.\n",d_lr);
	sasa_lee_richards(sasa,protein_xyz(p),r,protein_n(p),d_lr);
	break;
    default:
	fprintf(stderr,"Error: no SASA algorithm specified.\n");
	return 0;
    }

    printf("\n# SASA values in Å^2 for polar and non-polar atoms\n");
    sasa_per_atomclass(stdout,oons_classes(),p,sasa);

    printf("\n# SASA values in Å^2 for atom classes\n");
    sasa_per_atomclass(stdout,oons_types(),p,sasa);
    
    printf("\n# SASA values in Å^2 for residue types\n");
    sasa_per_atomclass(stdout,atomclassifier_residue(),p,sasa);
    
    printf("\n");
    
    protein_free(p);
    free(sasa);
    free(r);
    fclose(input);
}
