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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "freesasa.h"
#include "srp.h"

#if STDC_HEADERS
extern int getopt(int, char * const *, const char *);
extern int optind;
#endif

#ifdef PACKAGE_VERSION
const char *version = PACKAGE_VERSION;
#else
const char* version = "unknown";
#endif

char *program_name;

typedef struct  {
    freesasa *s;
    FILE *B;
    FILE *per_residue_type_file;
    FILE *per_residue_file;
    int per_residue_type;
    int per_residue;
    int printlog;
} settings_t;

static settings_t def_settings = {
    .s = NULL,
    .B = NULL,
    .per_residue_type_file = NULL,
    .per_residue_file = NULL,
    .per_residue_type = 0,
    .per_residue = 0,
    .printlog = 1,
};

static struct option long_options[] = {
    {"lee-richards", no_argument, 0, 'L'},
    {"shrake-rupley", no_argument, 0, 'S'},
    {"probe-radius", required_argument, 0, 'p'},
    {"lr-slice", required_argument, 0, 'd'},
    {"sr-points", required_argument, 0, 'n'},
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'v'},
    {"no-log", no_argument, 0, 'l'},
    {"no-warnings", no_argument, 0, 'w'},
    {"n-threads",required_argument, 0, 't'},
    {"hetatm",no_argument, 0, 'H'},
    {"sasa-per-residue-type",optional_argument,0,'r'},
    {"sasa-per-residue-sequence",optional_argument,0,'R'},
    {"print-as-B-values",required_argument,0,'B'}
};

void help() {
    fprintf(stderr,"\nUsage: %s [-hvlwLHSRrn:d:t:p:B:] pdb-file(s)\n",
            program_name);
    fprintf(stderr,
            "\nOptions are:\n"
            "  -h (--help)           Print this message\n"
            "  -v (--version)        Print version of the program\n"
            "  -l (--no-log)         Don't print log message\n"
            "  -w (--no-warnings)    Don't print warnings\n"
            "  -S (--shrake-rupley)  Use Shrake & Rupley algorithm [default]\n"
            "  -L (--lee-richards)   Use Lee & Richards algorithm\n"
            "  -H (--hetatm)         Include HETATM entries from input\n");
    
    fprintf(stderr,
            "  -p <value>  --probe-radius=<value>\n"
            "      Probe radius. Default value is %4.2f Å.\n"
            "  -d <value>  --lr-slice=<value>\n"
            "      Slice spacing in Lee & Richards algorithm. "
            "Default value is %4.2f Å\n"
            ,FREESASA_DEF_PROBE_RADIUS,FREESASA_DEF_LR_D);
    fprintf(stderr,
            "  -n <value>  --sr-points=<value>\n"
            "      Set number of test points in Shrake & Rupley algorithm\n"
            "      Default is %d, allowed values are:\n        ", FREESASA_DEF_SR_N);
    freesasa_srp_print_n_opt(stderr);
#ifdef HAVE_LIBPTHREAD
    fprintf(stderr,
            "  -t <value>  --n-threads=<value>\n"
            "      Set number of threads to use in calculation (for multicore CPUs)\n");
#endif
    fprintf(stderr,
            "  -r  --sasa-per-residue-type[=<output-file>]\n"
            "      Print SASA for each residue type. "
            "Writes to stdout if no output is specified.\n"
            "  -R  --sasa-per-residue-sequence[=<output-file>]\n"
            "      Print SASA for each residue, sequentially. "
            "Writes to stdout if no output is specified.\n");
    fprintf(stderr,
            "  -B  <output-file>  --print-as-B-values=<output-file>\n"
            "      Print PDB file with SASA for each atom as B-factors, "
            "in specified output file.\n");
    fprintf(stderr,
            "\nIf no pdb-file is specified STDIN is used for input.\n\n");
}

void short_help() {
    fprintf(stderr,"Run '%s -h' for usage instructions.\n",
            program_name);
}

void run_analysis(FILE *input, const char *name, const settings_t *settings) {
    freesasa *s = freesasa_new();
    double tmp;
    freesasa_copy_param(s,settings->s);
    freesasa_set_proteinname(s,name);
    if (freesasa_calc_pdb(s,input) == FREESASA_FAIL) {
        fprintf(stderr,"%s: error: Invalid input. Aborting.\n",
                program_name);
        exit(EXIT_FAILURE);
    }
    if (settings->printlog) {
        freesasa_log(s,stdout);
        printf("\nTotal:  %9.2f A2\nPolar:  %9.2f A2\nApolar: %9.2f A2\n",
               freesasa_area_total(s), freesasa_area_class(s,FREESASA_POLAR),
               freesasa_area_class(s, FREESASA_APOLAR));
    }
    if ((tmp = freesasa_area_class(s, FREESASA_NUCLEICACID)) > 0)
        printf("Nucleic: %9.2f A2\n",tmp);
    if ((tmp = freesasa_area_class(s, FREESASA_CLASS_UNKNOWN)) > 0)
        printf("Unknown: %9.2f A2\n",tmp);
    if (settings->per_residue_type) {
        FILE *f;
        if ((f = settings->per_residue_type_file) != NULL) {
            freesasa_per_residue_type(s,f);
            fclose(f);
        }
        else {
            freesasa_per_residue_type(s,stdout);
        }
    }
    if (settings->per_residue) {
        FILE *f;
        if ((f = settings->per_residue_file) != NULL) {
            freesasa_per_residue(s,f);
            fclose(f);
        } else {
            freesasa_per_residue(s,stdout);
        }
    }
    if (settings->B) freesasa_write_pdb(s,settings->B);
    freesasa_free(s);
}

int main (int argc, char **argv) {
    int alg_set = 0;
    FILE *input = NULL;
    settings_t settings = def_settings;
    settings.s = freesasa_new();
    extern char *optarg;
    char opt;
    int n_opt = 'z';
    char opt_set[n_opt];
    memset(opt_set,0,n_opt);
#ifdef PACKAGE_NAME
    program_name = PACKAGE_NAME;
#else
    program_name = argv[0];
#endif
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "hvlwLSHRrn:d:t:p:B:",
                              long_options, &option_index)) != -1) {
        errno = 0;
        opt_set[opt] = 1;
        switch(opt) {
        case 'h':
            help();
            exit(EXIT_SUCCESS);
        case 'v':
            printf("%s\n",version);
            exit(EXIT_SUCCESS);
        case 'l':
            settings.printlog = 0;
            break;
        case 'w':
            freesasa_set_verbosity(FREESASA_V_NOWARNINGS);
            break;
        case 'B':
            settings.B = fopen(optarg, "w");
            if (settings.B == NULL) {
                fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                        program_name,optarg,strerror(errno));
                short_help();
                exit(EXIT_FAILURE);
            }
            break;
        case 'n':
            freesasa_set_sr_points(settings.s,atoi(optarg));
            break;
        case 'S':
            freesasa_set_algorithm(settings.s,FREESASA_SHRAKE_RUPLEY);
            ++alg_set;
            break;
        case 'L':
            freesasa_set_algorithm(settings.s,FREESASA_LEE_RICHARDS);
            ++alg_set;
            break;
        case 'd':
            freesasa_set_lr_delta(settings.s,atof(optarg));
            break;
        case 'p':
            freesasa_set_probe_radius(settings.s,atof(optarg));
            break;
        case 'H':
            freesasa_include_hetatm(settings.s,1);
            break;
        case 'r':
            settings.per_residue_type = 1;
            if (optarg) {
                settings.per_residue_type_file = fopen(optarg,"w");
                if (settings.per_residue_type_file == NULL) {
                    fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                            program_name,optarg,strerror(errno));
                    short_help();
                    exit(EXIT_FAILURE);
                }
            }
            break;
        case 'R':
            settings.per_residue = 1;
            if (optarg) {
                settings.per_residue_file = fopen(optarg,"w");
                if (settings.per_residue_type_file == NULL) {
                    fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                            program_name,optarg,strerror(errno));
                    short_help();
                    exit(EXIT_FAILURE);
                }
            }
            break;
        case 't':
#if HAVE_LIBPTHREAD
            freesasa_set_nthreads(settings.s,atoi(optarg));
#else
            fprintf(stderr, "%s: warning: option 't' only defined if program"
                    " compiled with thread support.\n",
                    program_name);
#endif
            break;
        default:
            fprintf(stderr, "%s: warning: unknown option '%c' (will be ignored)\n",
                    program_name,opt);
            break;
        }
    }
    if (alg_set > 1) {
        fprintf(stderr, "%s: error: multiple algorithms specified.\n",
                program_name);
        exit(EXIT_FAILURE);
    }
    if ((opt_set['L'] && opt_set['n']) ||
        (!opt_set['L'] && opt_set['d']) ) {
        fprintf(stderr, "%s: warning: The program was given parameters "
                "not compatible with the selected algorithm. These will be ignored.\n",
                program_name);
    }

    if (argc > optind) {
        for (int i = optind; i < argc; ++i) {
            errno = 0;
            input = fopen(argv[i],"r");
            if (input != NULL) {
                run_analysis(input,argv[i],&settings);
                fclose(input);
            } else {
                fprintf(stderr, "%s: error opening file '%s'; %s\n",
                        program_name,argv[i],strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
    } else {
        if (!isatty(STDIN_FILENO)) {
            run_analysis(stdin,"stdin",&settings);
        } else {
            fprintf(stderr,"%s: no input.\n",
                    program_name);
            short_help();
        }
    }

    freesasa_free(settings.s);
    if (settings.B) fclose(settings.B);
    return EXIT_SUCCESS;
}
