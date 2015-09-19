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
#include <stdarg.h>
#if HAVE_CONFIG_H
#  include <config.h>
#endif

#include "freesasa.h"
#include "srp.h"

#if STDC_HEADERS
extern int getopt(int, char * const *, const char *);
extern int optind;
extern char *optarg;
#endif

#ifdef PACKAGE_VERSION
const char *version = PACKAGE_VERSION;
#else
const char* version = "unknown";
#endif

char *program_name;
const char* options_string;

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
int n_chain_groups = 0;
const char** chain_groups;

void
help(void)
{
    fprintf(stderr,"\nUsage: %s [%s] pdb-file(s)\n",
            program_name,options_string);
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
            "                        Number of threads to use in calculation. [default %d]\n",
            FREESASA_DEF_NUMBER_THREADS);
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
            "  -M (--separate-models) Calculate SASA for each MODEL separately.\n"
            "  -g <chains> --chain-groups <chains>\n"
            "                        Select chain or group of chains to treat separately.\n"
            "                        Several groups can be concatenated by '+', or by repetition.\n"
            "                        Examples:\n"
            "                                    '-g A', '-g AB', -g 'A+B', '-g A -g B', '-g AB+CD', etc\n");
    fprintf(stderr,"\nOutput options:\n"
            "  -l (--no-log)         Don't print log message (useful with -r -R and -B)\n"
            "  -w (--no-warnings)    Don't print warnings (will still print warnings due to invalid command\n"
            "                        line options)\n\n"
            "  -r  --per-residue-type --per-residue-type-file <output-file>\n"
            "  -R  --per-sequence --per-sequence-file <output-file>\n"
            "                        Print SASA for each residue, either grouped by type or sequentially.\n"
            "                        Use the -file variant to specify an output file.\n\n"
            "  -B  --print-as-B-values --B-value-file <output-file>\n"
            "                        Print PDB file where the temperature factor of each atom has\n"
            "                        been replaced by its SASA, and the occupancy number by the atomic\n"
            "                        radius. Use the -file variant to specify an output file.\n"
            "                        This option might at moment give confusing output when used in conjuction\n"
            "                        with the options -C and -g.\n");
    fprintf(stderr,
            "\nIf no pdb-file is specified STDIN is used for input.\n\n"
            "To calculate SASA of one or several PDB file using default parameters simply type:\n\n"
            "   '%s pdb-file(s)'     or    '%s < pdb-file'\n\n",
            program_name,program_name);
}

void
short_help(void)
{
    fprintf(stderr,"Run '%s -h' for usage instructions.\n",
            program_name);
}

void
abort_msg(const char *format,
          ...)
{
    va_list arg;
    va_start(arg, format);
    fprintf(stderr, "%s: error: ", program_name);
    vfprintf(stderr, format, arg);
    va_end(arg);
    fputc('\n', stderr);
    short_help();
    fputc('\n', stderr);
    fflush(stderr);
    exit(EXIT_FAILURE);
}


void
run_analysis(FILE *input,
const char *name) 
{
    double *radii;
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
    if (structures == NULL) abort_msg("Invalid input. Aborting.\n");

    if (n_chain_groups > 0) {
        int n2 = n;
        if (several_structures == 0) {
            structures = malloc(sizeof(freesasa_structure*));
            structures[0] = single_structure[0];
            several_structures = 1;
        }
        for (int i = 0; i < n_chain_groups; ++i) {
            for (int j = 0; j < n; ++j) {
                freesasa_structure* tmp = freesasa_structure_get_chains(structures[j],chain_groups[i]);
                if (tmp != NULL) {
                    ++n2;
                    structures = realloc(structures,sizeof(freesasa_structure*)*n2);
                    structures[n2-1] = tmp;
                } else {
                    abort_msg("Chain(s) '%s' not found.\n",chain_groups[i]);
                }
            }
        }
        n = n2;
    }
    for (int i = 0; i < n; ++i) {
        if (structures[i] == NULL) abort_msg("Invalid input.\n");
        radii = freesasa_structure_radius(structures[i],classifier);
        if (radii == NULL)         abort_msg("Can't calculate atomic radii.\n");
        result = freesasa_calc_structure(structures[i],radii,&parameters);
        if (result == NULL)       abort_msg("Can't calculate SASA.\n");
        classes = freesasa_result_classify(result,structures[i],classifier);
        if (classes == NULL)      abort_msg("Can't determine atom classes. Aborting.\n");
        if (printlog) {
            char name_i[name_len+10];
            strcpy(name_i,name);
            if (several_structures) {
                printf("\n");
                if (structure_options & FREESASA_SEPARATE_MODELS) 
                    sprintf(name_i+strlen(name_i),":%d",freesasa_structure_model(structures[i]));
                sprintf(name_i+strlen(name_i),":%s",freesasa_structure_chain_labels(structures[i]));
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

FILE*
fopen_werr(const char* filename,
           const char* mode) 
{
    errno = 0;
    FILE *f = fopen(filename,mode);
    if (f == NULL) {
        fprintf(stderr,"%s: error: could not open file '%s'; %s\n",
                program_name,optarg,strerror(errno));
        short_help();
        exit(EXIT_FAILURE);
    }
    return f;
}

void add_chain_groups(const char* cmd) 
{
    char *str = strdup(cmd);
    const char *token = strtok(str,"+");
    while (token) {
        ++n_chain_groups;
        chain_groups = realloc(chain_groups,sizeof(char*)*n_chain_groups);
        chain_groups[n_chain_groups-1] = strdup(token);
        token = strtok(0,",");
    }
    free(str);
}

int main (int argc, char **argv) 
{
    int alg_set = 0;
    FILE *input = NULL;
    char opt;
    int n_opt = 'z'+1;
    char opt_set[n_opt];
    int option_index = 0;
    static int option_flag;
    enum {B_FILE,RES_FILE,SEQ_FILE};
    parameters = freesasa_default_parameters;
    memset(opt_set,0,n_opt);
    // errors from this file will be prepended with freesasa, library errors with FreeSASA
    program_name = "freesasa";
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
        {"per-residue-type",no_argument,0,'r'},
        {"per-sequence",no_argument,0,'R'},
        {"print-as-B-values",no_argument,0,'B'},
        {"chain-groups",required_argument,0,'g'},
        {"per-residue-type-file",required_argument,&option_flag,RES_FILE},
        {"per-sequence-file",required_argument,&option_flag,SEQ_FILE},
        {"B-value-file",required_argument,&option_flag,B_FILE}
    };
    options_string = ":hvlwLSHYCMmBrRc:n:d:t:p:g:";
    while ((opt = getopt_long(argc, argv, options_string,
                              long_options, &option_index)) != -1) {
        opt_set[(int)opt] = 1;
        errno = 0;
        // Assume arguments starting with dash are actually missing arguments
        if (optarg != NULL && optarg[0] == '-') {
            if (option_index > 0) abort_msg("Missing argument? Value '%s' cannot be argument to '--%s'.\n",
                                            program_name,optarg,long_options[option_index].name);
            else abort_msg("Missing argument? Value '%s' cannot be argument to '-%c'.\n",
                           optarg,opt);
        }
        switch(opt) {
        case 0:
            switch(long_options[option_index].val) {
            case RES_FILE:
                per_residue_type = 1;
                per_residue_type_file = fopen_werr(optarg,"w");
                break;
            case SEQ_FILE:
                per_residue = 1;
                per_residue_file = fopen_werr(optarg,"w");
                break;
            case B_FILE:
                printpdb = 1;
                output_pdb = fopen_werr(optarg,"w");
                break;
            default:
                abort(); // what does this even mean?
            }
            break;
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
            FILE *f = fopen_werr(optarg,"r");
            classifier = freesasa_classifier_from_file(f);
            fclose(f);
            if (classifier == NULL) abort_msg("Can't read file '%s'.\n",optarg);
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
            if (parameters.probe_radius < 0) abort_msg("error: probe radius must be 0 or larger.\n");
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
            if (per_residue_type_file == NULL) per_residue_type_file = stdout;
            break;
        case 'R':
            per_residue = 1;
            if (per_residue_file == NULL) per_residue_file = stdout;
            break;
        case 'b':
        case 'B':
            printpdb = 1;
            if (output_pdb == NULL) output_pdb = stdout;
            break;
        case 'g':
            add_chain_groups(optarg);
            break;
        case 't':
#if HAVE_LIBPTHREAD
            parameters.n_threads = atoi(optarg);
            if (parameters.n_threads < 1) abort_msg("Number of threads must be 1 or larger.\n");
#else
            abort_msg("Option '-t' only defined if program compiled with thread support.\n");
#endif
            break;
        case ':':
            abort_msg("Option '-%c' missing argument.\n",optopt);
        case '?':
        default:
            fprintf(stderr, "%s: warning: Unknown option '-%c' (will be ignored)\n",
                    program_name,opt);
            break;
        }
    }
    if (alg_set > 1) abort_msg("Multiple algorithms specified.\n");
    if ((opt_set['L'] && opt_set['n']) ||
        (!opt_set['L'] && opt_set['d']) ) {
        abort_msg("The program was given parameters not compatible with the selected "
                  "algorithm. These will be ignored.\n");
    }
    if (opt_set['m'] && opt_set['M']) abort_msg("The options -m and -M can't be combined.\n");
    if (opt_set['g'] && opt_set['C']) abort_msg("The options -g and -C can't be combined.\n");
    if (printlog) printf("## %s %s ##\n",program_name,version);
    if (argc > optind) {
        for (int i = optind; i < argc; ++i) {
            errno = 0;
            input = fopen(argv[i],"r");
            if (input != NULL) {
                run_analysis(input,argv[i]);
                fclose(input);
            } else {
                abort_msg("Opening file '%s'; %s\n",argv[i],strerror(errno));
            }
        }
    } else {
        if (!isatty(STDIN_FILENO)) run_analysis(stdin,"stdin");
        else abort_msg("No input.\n", program_name);
    }

    if (classifier) freesasa_classifier_free(classifier);
    if (output_pdb) fclose(output_pdb);
    if (per_residue_type_file) fclose(per_residue_type_file);
    if (per_residue_file) fclose(per_residue_file);
    return EXIT_SUCCESS;
}
