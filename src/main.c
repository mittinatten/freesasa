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

freesasa_parameters parameters;
freesasa_classifier *classifier = NULL;
FILE *output_pdb = NULL;
FILE *per_residue_type_file = NULL;
FILE *per_residue_file = NULL;
int per_residue_type = 0;
int per_residue = 0;
int printlog = 1;
int structure_options = 0;
int printpdb = 0;

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
    {"hydrogen",no_argument, 0, 'Y'},
    {"separate-chains",no_argument, 0, 'C'},
    {"separate-models",no_argument, 0, 'M'},
    {"join-models",no_argument, 0, 'm'},
    {"config-file",required_argument,0,'c'},
    {"sasa-per-residue-type",optional_argument,0,'r'},
    {"sasa-per-residue-sequence",optional_argument,0,'R'},
    {"print-as-B-values",optional_argument,0,'B'}
};

void help() {
    fprintf(stderr,"\nUsage: %s [-hvlwLHYCMmSR::r::B::c:n:d:t:p:] pdb-file(s)\n",
            program_name);
    fprintf(stderr,"\n"
            "  -h (--help)           Print this message\n"
            "  -v (--version)        Print version of the program\n");
    fprintf(stderr,"\nSASA calculation parameters:\n"
            "  -S (--shrake-rupley)  Use Shrake & Rupley algorithm [default]\n"
            "  -L (--lee-richards)   Use Lee & Richards algorithm\n");
    fprintf(stderr,
            "  -p <value>  --probe-radius=<value>\n"
            "                        Probe radius [default %4.2f Å]\n"
            "  -d <value>  --lr-slice=<value>\n"
            "                        Slice spacing in Lee & Richards algorithm \n"
            "                        [default %4.2f Å].\n",
            FREESASA_DEF_PROBE_RADIUS,FREESASA_DEF_LR_D);
    fprintf(stderr,
            "  -n <value>  --sr-points=<value>\n"
            "                        Number of test points in Shrake & Rupley algorithm.\n"
            "                        Default is %d, allowed values are:\n"
            "                          ", FREESASA_DEF_SR_N);
    freesasa_srp_print_n_opt(stderr);
#ifdef HAVE_LIBPTHREAD
    fprintf(stderr,
            "  -t <value>  --n-threads=<value>\n"
            "                        Number of threads to use in calculation.\n");
#endif
    fprintf(stderr,
            "  -c <file> (--config-file=<file>)\n"
            "                        Use atomic radii and classes provided in file\n");
    fprintf(stderr,
            "\nInput PDB:              (default: ignore HETATM and hydrogens, include all\n"
            "                         chains of the first MODEL)\n"
            "  -H (--hetatm)         Include HETATM entries from input\n"
            "  -Y (--hydrogen)       Include hydrogen atoms. Only makes sense in concjunction\n"
            "                        with -c option. Default H radius is 0 Å.\n"
            "  -m (--join-models)    Join all MODELs in input into one big structure.\n"
            "  -C (--separate-chains) Calculate SASA for each chain separately.\n"
            "  -M (--separate-models) Calculate SASA for each MODEL separately.\n");
    fprintf(stderr,"\nOutput options:\n"
            "  -l (--no-log)         Don't print log message (useful with -r -R and -B)\n"
            "  -w (--no-warnings)    Don't print warnings\n"
            "  -r  --sasa-per-residue-type[=<output-file>]\n"
            "  -R  --sasa-per-residue-sequence[=<output-file>]\n"
            "                        Print SASA for each residue, either grouped by type or sequentially.\n"
            "                        Writes to STDOUT if no output is specified.\n"
            "  -B  --print-as-B-values[=<output-file>]\n"
            "                        Print PDB file with SASA for each atom as B-factors.\n"
            "                        Write to STDOUT if no output is specified.\n");
    fprintf(stderr,
            "\nIf no pdb-file is specified STDIN is used for input.\n\n");
}

void short_help() {
    fprintf(stderr,"Run '%s -h' for usage instructions.\n",
            program_name);
}

void run_analysis(FILE *input, const char *name) {
    double tmp, *radii;
    int several_structures = 0, name_len = strlen(name);
    freesasa_result *result;
    freesasa_strvp *classes = NULL;
    freesasa_structure *single_structure[1];
    freesasa_structure **structures;
    int n = 0;
    if ((structure_options & FREESASA_SEPARATE_CHAINS) ||
        (structure_options & FREESASA_SEPARATE_MODELS)) {
        structures = freesasa_structure_array(input,&n,structure_options);
        several_structures = 1;
    } else {
        single_structure[0] = freesasa_structure_from_pdb(input,structure_options);
        structures = single_structure;
        n = 1;
    }
    if (structures == NULL) {
        fprintf(stderr,"%s: error: Invalid input. Aborting.\n",
                program_name);
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < n; ++i) {
        if (structures[i] == NULL) {
            fprintf(stderr,"%s: error: Invalid input. Aborting.\n",
                    program_name);
            exit(EXIT_FAILURE);
        }
        radii = freesasa_structure_radius(structures[i],classifier);
        if (radii == NULL) {
            fprintf(stderr,"%s: error: Can't calculate atomic radii. Aborting.\n",
                    program_name);
            exit(EXIT_FAILURE);
        }
        if ((result = freesasa_calc_structure(structures[i],radii,&parameters)) == NULL) {
            fprintf(stderr,"%s: error: Can't calculate SASA. Aborting.\n",
                    program_name);
            exit(EXIT_FAILURE);
        }
        classes = freesasa_result_classify(result,structures[i],classifier);
        if (classes == NULL) {
            fprintf(stderr,"%s: warning: Can't determine atom classes. Aborting.\n",
                    program_name);
            exit(EXIT_FAILURE);
        }
        if (printlog) {
            char name_i[name_len+10];
            strcpy(name_i,name);
            if (several_structures) {
                printf("\n");
                if (structure_options & FREESASA_SEPARATE_MODELS) 
                    sprintf(name_i+strlen(name_i),":%d",freesasa_structure_model(structures[i]));
                if (structure_options & FREESASA_SEPARATE_CHAINS) 
                    sprintf(name_i+strlen(name_i),":%c",freesasa_structure_atom_chain(structures[i],0));
            }
            freesasa_log(stdout,result,name_i,&parameters,classes);
            if (several_structures) printf("\n");
        }
        if (per_residue_type) {
            if (several_structures) fprintf(per_residue_type_file,"\n## Structure %d\n",i);
            freesasa_per_residue_type(per_residue_type_file,result,structures[i]);
        }
        if (per_residue) {
            if (several_structures) fprintf(per_residue_file,"\n## Structure %d\n",i);
            freesasa_per_residue(per_residue_file,result,structures[i]);
        }
        if (printpdb) {
            freesasa_write_pdb(output_pdb,result,structures[i],radii);
        }
        freesasa_result_free(result);
        freesasa_strvp_free(classes);
        freesasa_structure_free(structures[i]);
        free(radii);
    }
    if (structures != single_structure) free(structures);
}

int main (int argc, char **argv) {
    int alg_set = 0;
    FILE *input = NULL;
    extern char *optarg;
    parameters = freesasa_default_parameters;
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
    while ((opt = getopt_long(argc, argv, "hvlwLSHYCMmR::r::B::c:n:d:t:p:",
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
            printlog = 0;
            break;
        case 'w':
            freesasa_set_verbosity(FREESASA_V_NOWARNINGS);
            break;
        case 'c': {
            FILE *f = fopen(optarg,"r");
            if (f == NULL) {
                fprintf(stderr,"%s: error: could not open file '%s'; %s\nAborting.",
                        program_name,optarg,strerror(errno));
                short_help();
                exit(EXIT_FAILURE);
            } else {
                classifier = freesasa_classifier_from_file(f);
                if (classifier == NULL) {
                    fclose(f);
                    fprintf(stderr,"%s: error: Can't read file '%s'. Aborting.\n",
                            program_name,optarg);
                    exit(EXIT_FAILURE);
                }
            }
            break;
        }
        case 'n':
            parameters.shrake_rupley_n_points = atoi(optarg);
            break;
        case 'S':
            parameters.alg = FREESASA_SHRAKE_RUPLEY;
            ++alg_set;
            break;
        case 'L':
            parameters.alg = FREESASA_LEE_RICHARDS;
            ++alg_set;
            break;
        case 'd':
            parameters.lee_richards_delta = atof(optarg);
            break;
        case 'p':
            parameters.probe_radius = atof(optarg);
            if (parameters.probe_radius < 0) {
                fprintf(stderr, "%s: error: probe radius must be 0 or larger.\n",
                        program_name);
                exit(EXIT_FAILURE);
            }
            break;
        case 'H':
            structure_options |= FREESASA_INCLUDE_HETATM;
            break;
        case 'Y':
            structure_options |= FREESASA_INCLUDE_HYDROGEN;
            break;
        case 'M':
            structure_options |= FREESASA_SEPARATE_MODELS;
            break;
        case 'm':
            structure_options |= FREESASA_JOIN_MODELS;
            break;
        case 'C':
            structure_options |= FREESASA_SEPARATE_CHAINS;
            break;
        case 'r':
            per_residue_type = 1;
            if (optarg) {
                per_residue_type_file = fopen(optarg,"w");
                if (per_residue_type_file == NULL) {
                    fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                            program_name,optarg,strerror(errno));
                    short_help();
                    exit(EXIT_FAILURE);
                }
            } else {
                per_residue_type_file = stdout;
            }
            break;
        case 'R':
            per_residue = 1;
            if (optarg) {
                per_residue_file = fopen(optarg,"w");
                if (per_residue_file == NULL) {
                    fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                            program_name,optarg,strerror(errno));
                    short_help();
                    exit(EXIT_FAILURE);
                }
            } else {
                per_residue_file = stdout;
            }
            break;
        case 'B':
            printpdb = 1;
            if (optarg) {
                output_pdb = fopen(optarg,"w");
                if (output_pdb == NULL) {
                    fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                            program_name,optarg,strerror(errno));
                    short_help();
                    exit(EXIT_FAILURE);
                }
            } else {
                output_pdb = stdout;
            }
            break;
        case 't':
#if HAVE_LIBPTHREAD
            parameters.n_threads = atoi(optarg);
            if (parameters.n_threads < 1) {
                fprintf(stderr, "%s: error: number of threads must be 1 or larger.\n",
                        program_name);
                exit(EXIT_FAILURE);
            }
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
    if (opt_set['m'] && opt_set['M']) {
        fprintf(stderr, "%s: error: The options -m and -M can't be combined.\n",
                program_name);
        exit(EXIT_FAILURE);
    }
    if (printlog) printf("## %s %s ##\n",program_name,version);
    if (argc > optind) {
        for (int i = optind; i < argc; ++i) {
            errno = 0;
            input = fopen(argv[i],"r");
            if (input != NULL) {
                run_analysis(input,argv[i]);
                fclose(input);
            } else {
                fprintf(stderr, "%s: error opening file '%s'; %s\n",
                        program_name,argv[i],strerror(errno));
                exit(EXIT_FAILURE);
            }
        }
    } else {
        if (!isatty(STDIN_FILENO)) {
            run_analysis(stdin,"stdin");
        } else {
            fprintf(stderr,"%s: no input.\n",
                    program_name);
            short_help();
            exit(EXIT_FAILURE);
        }
    }

    if (classifier) freesasa_classifier_free(classifier);
    if (output_pdb) fclose(output_pdb);
    if (per_residue_type_file) fclose(per_residue_type_file);
    if (per_residue_file) fclose(per_residue_file);
    return EXIT_SUCCESS;
}
