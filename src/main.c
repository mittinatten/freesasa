#if HAVE_CONFIG_H
#  include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <getopt.h>
#include <stdarg.h>

#include "freesasa.h"
#include "freesasa_internal.h"

#if STDC_HEADERS
extern int getopt(int, char * const *, const char *);
extern int optind, optopt;
extern char *optarg;
#endif

char *program_name = "freesasa";
const char* options_string;

// configuration
freesasa_parameters parameters;
const freesasa_classifier *classifier = NULL;
freesasa_classifier *classifier_from_file = NULL;
int structure_options = 0;

// output files
FILE *output = NULL;
FILE *errlog;

// flags
int printlog = 1;
int static_config = 0;
int output_depth = FREESASA_OUTPUT_CHAIN;
int no_rel = 0;

// chain groups
int n_chain_groups = 0;
char** chain_groups = NULL;

// selection commands
int n_select = 0;
char** select_cmd = NULL;

struct analysis_results {
    freesasa_node *tree;
    freesasa_structure **structures;
    int n_structures;
};

static void
help(void)
{
    fprintf(stderr, "\nUsage: %s [options] pdb-file(s)\n\n", program_name);
    fprintf(stderr, "GENERAL OPTIONS\n"
            "  -h (--help)           Print this message\n"
            "  -v (--version)        Print version of the program\n"
            "  --deprecated          Print deprecated options\n");
    fprintf(stderr, "\nPARAMETERS\n"
            "  -S (--shrake-rupley)  Use Shrake & Rupley algorithm\n"
            "  -L (--lee-richards)   Use Lee & Richards algorithm [default]\n");
    fprintf(stderr,
            "\n"
            "  -p <value>  (--probe-radius=<value>)\n"
            "                        Probe radius [default: %4.2f Å]\n"
            "\n"
            "  -n <value>  (--resolution=<value>)\n"
            "                        Either: \n"
            "                        - Number of test points in Shrake & Rupley algorithm,\n"
            "                          [default: %d] or\n"
            "                        - number of slices per atom in Lee & Richards algorithm.\n"
            "                          [default: %d]\n"
            "                        depending on which is selected.\n",
            FREESASA_DEF_PROBE_RADIUS, FREESASA_DEF_SR_N, FREESASA_DEF_LR_N);
    if (USE_THREADS) {
        fprintf(stderr,
                "\n  -t <value>  (--n-threads=<value>)\n"
                "                        Number of threads to use in calculation. [default %d]\n",
                FREESASA_DEF_NUMBER_THREADS);
    }
    fprintf(stderr,
            "\n  -O (--radius-from-occupancy)\n"
            "                        Read atomic radii from Occupancy field in the PDB input.\n"
            "\n  -c <file> (--config-file=<file>)\n"
            "                        Use atomic radii and classes provided in file, example\n"
            "                        configuration files can be found in the directory\n"
            "                        share/.\n"
            "\n  --radii=<protor|naccess>\n"
            "                        Use either ProtOr or NACCESS atomic radii, classes and\n"
            "                        RSA reference values. Cannot be used in conjunction\n"
            "                        with the option '-c'.\n"
            "                        Default value is 'protor'.\n");
    fprintf(stderr, "\nINPUT\n"
            "  -H (--hetatm)         Include HETATM entries from input.\n"
            "  -Y (--hydrogen)       Include hydrogen atoms (skipped by default). Default\n"
            "                        classifier emits warnings. Use with care. To get\n"
            "                        sensible results, one probably needs to redefine atomic\n"
            "                        radii with the -c option. Default H radius is 1.10 Å.\n"
            "  -m (--join-models)    Join all MODELs in input into one big structure.\n"
            "  -C (--separate-chains) Calculate SASA for each chain separately.\n"
            "  -M (--separate-models) Calculate SASA for each MODEL separately.\n"
            "\n"
            "  -g <chains> (--chain-groups=<chains>)\n"
            "                        Select chain or group of chains, several groups can be \n"
            "                        concatenated by '+', or by repeting the command. A\n"
            "                        separate SASA calculation will be performed for each\n"
            "                        group as though the other chains didn't exist.\n\n"
            "                        Examples:\n"
            "                            '-g A', '-g AB', -g 'A+B', '-g A -g B', '-g AB+CD'\n"
            "\n"
            "  --unknown=<guess|skip|halt>\n"
            "                        When an unknown atom is encountered FreeSASA can either\n"
            "                        'guess' its VdW radius, 'skip' the atom, or 'halt'.\n"
            "                        Default is 'guess'.\n");
    fprintf(stderr, "\nOUTPUT\n"
            "  -l (--no-log)         Don't print log message (useful with -r and -R)\n"
            "  -w (--no-warnings)    Don't print warnings (will still print warnings due to\n"
            "                        invalid command line options)\n"
            "\n"
            "  -o <file> (--output=<file>)\n"
            "  -e <file> (--error-file=<file>)\n"
            "                        Redirect output and/or errors and warnings to file.\n"
            "\n"
            "  -f <format> (--format <format>)\n"
            "                        Output format, options are:\n"
            "                          - log  Default output, plain text\n"
            "                          - res  The SASA of each residue type on a separate line.\n"
            "                          - seq  The SASA of each residue separately.\n"
            "                          - pdb  PDB file where the temperature factor of each atom\n"
            "                                 is replaced by its SASA, and the occupancy number\n"
            "                                 by the atomic radius.\n"
            "                          - rsa  Results in RSA format.\n"
            "                          - json Results in JSON format.\n"
            "                          - xml  Results in XML format.\n"
            "                        If the option is repeated, only the last value will be used.\n"
            "\n");
    if (USE_JSON || USE_XML) {
        fprintf(stderr,
                "  --output-depth=<structure|chain|residue|atom>\n"
                "                        Level of detail in JSON and/or XML output [default: chain]\n"
                "\n");
    }
    fprintf(stderr,
            "\n"
            "  --select <command>    Select atoms using Pymol select syntax.\n"
            "                        The option can be repeated to define several selections.\n\n"
            "                        Examples:\n"
            "                            'AR, resn ala+arg', 'chain_A, chain A'\n"
            "                        AR and chain_A are just the names of the selections,\n"
            "                        which will be reused in output. See documentation for\n"
            "                        full syntax specification.\n");
    fprintf(stderr,
            "\nIf no pdb-file is specified STDIN is used for input.\n\n"
            "To calculate SASA of one or several PDB file using default parameters simply\ntype:\n\n"
            "   '%s pdb-file(s)'     or    '%s < pdb-file'\n\n",
            program_name,program_name);
}

static void
deprecated(void)
{
    fprintf(stderr,
            "These options are still supported but will disappear in later versions of FreeSASA.\n"
            "Use --format instead\n\n"
            "  --rsa                 Equivalent to --format=rsa\n"
            "  -B  (--print-as-B-values)\n"
            "                        Equivalent to --format=pdb\n"
            "  -r  (--foreach-residue-type)\n"
            "                        Equivalent to --format=res\n"
            "  -R  (--foreach-residue)\n"
            "                        Equivalent to --format=seq.\n"

            "\n");
}

static void
short_help(void)
{
    fprintf(stderr, "Run '%s -h' for usage instructions.\n",
            program_name);
}

static void
version(void)
{
    printf("%s\n", freesasa_string);
    printf("License: MIT <http://opensource.org/licenses/MIT>\n");
    printf("If you use this program for research, please cite:\n");
    printf("  Simon Mitternacht (2016) FreeSASA: An open source C\n"
           "  library for solvent accessible surface area calculations.\n"
           "  F1000Research 5:189.\n");
}

static void
release_resources(void)
{
    if (classifier_from_file) freesasa_classifier_free(classifier_from_file);
    if (errlog) fclose(errlog);
    if (chain_groups) {
        for (int i = 0; i < n_chain_groups; ++i) {
            free(chain_groups[i]);
        }
    }
    if (select_cmd) {
        for (int i = 0; i < n_select; ++i) {
            free(select_cmd[i]);
        }
    }
}

static void
abort_msg(const char *format,
          ...)
{
    va_list arg;
    va_start(arg, format);
    fprintf(stderr, "%s: error: ", program_name);
    vfprintf(stderr, format, arg);
    va_end(arg);
    fputc('\n', stderr);
    fputc('\n', stderr);
    short_help();
    fputc('\n', stderr);
    fflush(stderr);
    release_resources();
    exit(EXIT_FAILURE);
}

static freesasa_structure **
get_structures(FILE *input, int *n)
{
   freesasa_structure **structures = NULL;

   *n = 0;
   if ((structure_options & FREESASA_SEPARATE_CHAINS) ||
       (structure_options & FREESASA_SEPARATE_MODELS)) {
       structures = freesasa_structure_array(input, n, classifier, structure_options);
       if (structures == NULL) abort_msg("Invalid input.");
       for (int i = 0; i < *n; ++i) {
           if (structures[i] == NULL) abort_msg("Invalid input.");
       }
   } else {
       structures = malloc(sizeof(freesasa_structure*));
       if (structures == NULL) {
           abort_msg("Out of memory.");
       }
       *n = 1;
       structures[0] = freesasa_structure_from_pdb(input, classifier, structure_options);
       if (structures[0] == NULL) {
           abort_msg("Invalid input.");
       }
   }

   // get chain-groups (if requested)
   if (n_chain_groups > 0) {
       int n2 = *n;
       for (int i = 0; i < n_chain_groups; ++i) {
           for (int j = 0; j < *n; ++j) {
               freesasa_structure* tmp = freesasa_structure_get_chains(structures[j], chain_groups[i]);
               if (tmp != NULL) {
                   ++n2;
                   structures = realloc(structures, sizeof(freesasa_structure*)*n2);
                   if (structures == NULL) abort_msg("Out of memory.");
                   structures[n2-1] = tmp;
               } else {
                   abort_msg("Chain(s) '%s' not found.", chain_groups[i]);
               }
           }
       }
       *n = n2;
   }

   return structures;
}

static freesasa_node *
run_analysis(FILE *input,
             const char *name)
{
    int name_len = strlen(name);
    freesasa_structure **structures = NULL;
    freesasa_node *tree = freesasa_tree_new();
    int n = 0;

    if (tree == NULL) abort_msg("Failed to initialize result-tree.");

    // read PDB file
    structures = get_structures(input, &n);
    if (n == 0) abort_msg("Invalid input.");

    // perform calculation on each structure and output results
    for (int i = 0; i < n; ++i) {
        char name_i[name_len+10];
        freesasa_node *tmp_tree;

        strcpy(name_i,name);
        if (n > 1 && (structure_options & FREESASA_SEPARATE_MODELS))
            sprintf(name_i+strlen(name_i), ":%d", freesasa_structure_model(structures[i]));

        tmp_tree = freesasa_calc_tree(structures[i], &parameters, name_i);
        if (tmp_tree == NULL) abort_msg("Can't calculate SASA.");

        freesasa_node *structure_node =
            freesasa_node_children(freesasa_node_children(tmp_tree));
        const freesasa_result *result = freesasa_node_structure_result(structure_node);

        // Calculate selections for each structure
        if (n_select > 0) {
            for (int c = 0; c < n_select; ++c) {
                freesasa_selection *sel = freesasa_selection_new(select_cmd[c], structures[i], result);
                if (sel != NULL) {
                    freesasa_node_structure_add_selection(structure_node, sel);
                } else {
                    abort_msg("Illegal selection");
                }
                freesasa_selection_free(sel);
            }
        }

        if (freesasa_tree_join(tree, &tmp_tree) != FREESASA_SUCCESS) {
            abort_msg("Failed joining result-trees");
        }

        freesasa_structure_free(structures[i]);
    }

    free(structures);

    return tree;
}

static FILE*
fopen_werr(const char* filename,
           const char* mode) 
{
    FILE *f = fopen(filename, mode);
    if (f == NULL) {
        abort_msg("could not open file '%s'; %s",
                  filename, strerror(errno));
    }
    return f;
}

static void
add_chain_groups(const char* cmd) 
{
    int err = 0;
    char *str;
    const char *token;

    // check that string is valid
    for (size_t i = 0; i < strlen(cmd); ++i) {
        char a = cmd[i];
        if (a != '+' && 
            !(a >= 'a' && a <= 'z') && !(a >= 'A' && a <= 'Z') &&
            !(a >= '0' && a <= '9')) {
            freesasa_fail("Character '%c' not valid chain ID in --chain-groups. "
                          "Valdig characters are [A-z0-9] and '+' as separator.",a);
            ++err;
        }
    }

    //extract chain groups
    if (err == 0) {
        str = strdup(cmd);
        token = strtok(str, "+");
        while (token) {
            ++n_chain_groups;
            chain_groups = realloc(chain_groups,sizeof(char*)*n_chain_groups);
            if (chain_groups == NULL) { mem_fail(); abort();}
            chain_groups[n_chain_groups-1] = strdup(token);
            token = strtok(0, "+");
        }
        free(str);
    } else {
        abort_msg("Aborting.");
    }
}

static void
add_select(const char* cmd) 
{
    ++n_select;
    select_cmd = realloc(select_cmd, sizeof(char*)*n_select);
    if (select_cmd == NULL) {
        mem_fail(); 
        abort(); 
    }
    select_cmd[n_select-1] = strdup(cmd);
    if (select_cmd[n_select-1] == NULL) {
        mem_fail(); 
        abort();
    }
}

static void
add_unknown_option(const char *optarg)
{
    if (strcmp(optarg, "skip") == 0) {
        structure_options |= FREESASA_SKIP_UNKNOWN;
        return;
    } 
    if(strcmp(optarg, "halt") == 0) {
        structure_options |= FREESASA_HALT_AT_UNKNOWN;
        return;
    } 
    if(strcmp(optarg, "guess") == 0) {
        return; //default
    }
    abort_msg("Unknown alternative to option --unknown: '%s'", optarg);
}

static int
parse_output_format(const char *optarg)
{
    if (strcmp(optarg, "log") == 0) {
        return FREESASA_LOG;
    }
    if (strcmp(optarg, "res") == 0) {
        return FREESASA_RES;
    }
    if (strcmp(optarg, "seq") == 0) {
        return FREESASA_SEQ;
    } 
    if (strcmp(optarg, "rsa") == 0) {
        return FREESASA_RSA;
    }
    if (strcmp(optarg, "json") == 0) {
        return FREESASA_JSON;
    }
    if (strcmp(optarg, "xml") == 0) {
        return FREESASA_XML;
    }
    if (strcmp(optarg, "pdb") == 0) {
        return FREESASA_PDB;
    }
    abort_msg("Unknown output format: '%s'", optarg);
    return FREESASA_FAIL; // to avoid compiler warnings
}

static int
parse_output_depth(const char *optarg) {
    if (strcmp("structure", optarg) == 0) {
        return FREESASA_OUTPUT_STRUCTURE;
    }
    if (strcmp("chain", optarg) == 0) {
        return FREESASA_OUTPUT_CHAIN;
    }
    if (strcmp("residue", optarg) == 0) {
        return FREESASA_OUTPUT_RESIDUE;
    }
    if (strcmp("atom", optarg) == 0) {
        return FREESASA_OUTPUT_ATOM;
    }
    abort_msg("Output depth '%s' not allowed, "
              "can only be 'structure', 'chain', 'residue' or 'atom'",
              optarg);
    return FREESASA_FAIL; // to avoid compiler warnings
}

static const freesasa_classifier *
parse_radii_option(const char *optarg)
{
    if (strcmp("naccess", optarg) == 0) {
        return &freesasa_naccess_classifier;
    } else if (strcmp("protor", optarg) == 0) {
        return &freesasa_protor_classifier;
    }
    abort_msg("Config '%s' not allowed, "
              "can only be 'protor' or 'naccess')", optarg);
    return NULL; // to avoid compiler warnings
}

int
main(int argc,
     char **argv) 
{
    int alg_set = 0;
    FILE *input = NULL;
    char opt;
    int n_opt = 'z'+1;
    char opt_set[n_opt];
    int option_index = 0;
    int option_flag;
    int output_format = FREESASA_LOG;
    freesasa_node *tree = freesasa_tree_new();
    enum {B_FILE, SELECT, UNKNOWN, RSA, JSON, XML,
          O_DEPTH, RADII, DEPRECATED};
    parameters = freesasa_default_parameters;
    memset(opt_set, 0, n_opt);
    if (tree == NULL) abort_msg("Error initializing calculation.");

    struct option long_options[] = {
        {"lee-richards",         no_argument,       0, 'L'},
        {"shrake-rupley",        no_argument,       0, 'S'},
        {"probe-radius",         required_argument, 0, 'p'},
        {"resolution",           required_argument, 0, 'n'},
        {"help",                 no_argument,       0, 'h'},
        {"version",              no_argument,       0, 'v'},
        {"no-log",               no_argument,       0, 'l'},
        {"no-warnings",          no_argument,       0, 'w'},
        {"n-threads",            required_argument, 0, 't'},
        {"hetatm",               no_argument,       0, 'H'},
        {"hydrogen",             no_argument,       0, 'Y'},
        {"separate-chains",      no_argument,       0, 'C'},
        {"separate-models",      no_argument,       0, 'M'},
        {"join-models",          no_argument,       0, 'm'},
        {"config-file",          required_argument, 0, 'c'},
        {"foreach-residue-type", no_argument,       0, 'r'},
        {"foreach-residue",      no_argument,       0, 'R'},
        {"print-as-B-values",    no_argument,       0, 'B'},
        {"chain-groups",         required_argument, 0, 'g'},
        {"error-file",           required_argument, 0, 'e'},
        {"output",               required_argument, 0, 'o'},
        {"format",               required_argument, 0, 'f'},
        {"radius-from-occupancy",no_argument,       0, 'O'},
        {"select",               required_argument, &option_flag, SELECT},
        {"unknown",              required_argument, &option_flag, UNKNOWN},
        {"rsa",                  no_argument,       &option_flag, RSA},
        {"output-depth",         required_argument, &option_flag, O_DEPTH},
        {"radii",                required_argument, &option_flag, RADII},
        {"deprecated",           no_argument,       &option_flag, DEPRECATED},
        {0,0,0,0}
    };
    options_string = ":hvlwLSHYOCMmBrRc:n:t:p:g:e:o:f:";
    while ((opt = getopt_long(argc, argv, options_string,
                              long_options, &option_index)) != -1) {
        opt_set[(int)opt] = 1;
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
            case SELECT:
                add_select(optarg);
                break;
            case UNKNOWN:
                add_unknown_option(optarg);
                break;
            case RSA:
                output_format = FREESASA_RSA;
                break;
            case O_DEPTH:
                output_depth = parse_output_depth(optarg);
                break;
            case RADII:
                static_config = 1;
                classifier = parse_radii_option(optarg);
                break;
            case DEPRECATED:
                deprecated();
                exit(EXIT_SUCCESS);
            default:
                abort(); // what does this even mean?
            }
            break;
        case 'h':
            help();
            exit(EXIT_SUCCESS);
        case 'v':
            version();
            exit(EXIT_SUCCESS);
        case 'e': 
            errlog = fopen_werr(optarg, "w");
            freesasa_set_err_out(errlog);
            break;
        case 'o':
            output = fopen_werr(optarg, "w");
            break;
        case 'f':
            output_format = parse_output_format(optarg);
            break;
        case 'l':
            printlog = 0;
            break;
        case 'w':
            freesasa_set_verbosity(FREESASA_V_NOWARNINGS);
            break;
        case 'c': {
            FILE *cf = fopen_werr(optarg, "r");
            classifier = classifier_from_file = freesasa_classifier_from_file(cf);
            if (classifier_from_file == NULL) abort_msg("Can't read file '%s'.", optarg);
            no_rel = 1;
            break;
        }
        case 'n':
            parameters.shrake_rupley_n_points = atoi(optarg);
            parameters.lee_richards_n_slices = atoi(optarg);
            if (parameters.shrake_rupley_n_points <= 0)
                abort_msg("Resolution needs to be at least 1 (20 recommended minum for S&R, 5 for L&R).");
            break;
        case 'S':
            parameters.alg = FREESASA_SHRAKE_RUPLEY;
            ++alg_set;
            break;
        case 'L':
            parameters.alg = FREESASA_LEE_RICHARDS;
            ++alg_set;
            break;
        case 'p':
            parameters.probe_radius = atof(optarg);
            if (parameters.probe_radius <= 0)
                abort_msg("Probe radius must be 0 or larger.");
            break;
        case 'H':
            structure_options |= FREESASA_INCLUDE_HETATM;
            break;
        case 'Y':
            structure_options |= FREESASA_INCLUDE_HYDROGEN;
            break;
        case 'O':
            structure_options |= FREESASA_RADIUS_FROM_OCCUPANCY;
            no_rel = 1;
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
            freesasa_warn("Option '-r' deprecated, use '-f res' or '--format=res' instead");
            output_format = FREESASA_RES;
            break;
        case 'R':
            freesasa_warn("Option '-R' deprecated, use '-f seq' or '--format=seq' instead");
            output_format = FREESASA_SEQ;
            break;
        case 'B':
            freesasa_warn("Option '-B' deprecated, use '-f pdb' or '--format=pdb' instead");
            output_format = FREESASA_PDB;
            break;
        case 'g':
            add_chain_groups(optarg);
            break;
        case 't':
            if (USE_THREADS) {
                parameters.n_threads = atoi(optarg);
                if (parameters.n_threads < 1) abort_msg("Number of threads must be 1 or larger.");
            } else {
                abort_msg("Option '-t' only defined if program compiled with thread support.");
            }
            break;
        case ':':
            abort_msg("Option '-%c' missing argument.", optopt);
            break;
        case '?':
        default:
            abort_msg("Unknown option '-%c'.", optopt);
            break;
        }
    }
    if (output == NULL) output = stdout;
    if (alg_set > 1) abort_msg("Multiple algorithms specified.");
    if (opt_set['m'] && opt_set['M']) abort_msg("The options -m and -M can't be combined.");
    if (opt_set['g'] && opt_set['C']) abort_msg("The options -g and -C can't be combined.");
    if (opt_set['c'] && static_config) abort_msg("The options -c and --radii cannot be combined");
    if (opt_set['O'] && static_config) abort_msg("The options -O and --radii cannot be combined");
    if (opt_set['c'] && opt_set['O']) abort_msg("The options -c and -O can't be combined");
    if (output_format == FREESASA_RSA && (opt_set['c'] || opt_set['O'])) {
        freesasa_warn("Will skip REL columns in RSA when custom atomic radii selected.");
    }
    if (output_format == FREESASA_RSA && (opt_set['C'] || opt_set['M']))
        abort_msg("The RSA format can not be used with the options -C or -M. "
                  "The format does not support several results in one file.");
    if (output_format == FREESASA_LOG) {
        fprintf(output, "## %s ##\n", freesasa_string);
    }
    if (argc > optind) {
        for (int i = optind; i < argc; ++i) {
            freesasa_node *tmp;
            input = fopen_werr(argv[i], "r");
            if (output_format == FREESASA_LOG) {
                freesasa_write_parameters(output, &parameters);
            }
            tmp = run_analysis(input, argv[i]);
            freesasa_tree_join(tree, &tmp);
            fclose(input);
        }
    } else {
        if (!isatty(STDIN_FILENO)) {
            freesasa_node *tmp;
            tmp = run_analysis(stdin, "stdin");
            freesasa_tree_join(tree, &tmp);
        }
        else abort_msg("No input.", program_name);
    }

    freesasa_tree_export(output, tree, output_format | output_depth | (no_rel ? FREESASA_OUTPUT_SKIP_REL : 0));
    freesasa_node_free(tree);

    release_resources();

    return EXIT_SUCCESS;
}
 
