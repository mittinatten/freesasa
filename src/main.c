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

#ifdef PACKAGE_VERSION
const char *version = PACKAGE_VERSION;
#else
const char* version = "unknown";
#endif

char *program_name;
const char* options_string;

// configuration
freesasa_parameters parameters;
const freesasa_classifier *classifier = NULL;
freesasa_classifier *classifier_from_file = NULL;
int structure_options = 0;

// output files
FILE *output_pdb = NULL;
FILE *per_residue_type_file = NULL;
FILE *per_residue_file = NULL;
FILE *rsa_file = NULL;
FILE *json_file = NULL;
FILE *xml_file = NULL;
FILE *output = NULL;
FILE *errlog;

// flags
int per_residue_type = 0;
int per_residue = 0;
int printlog = 1;
int printpdb = 0;
int printrsa = 0;
int printjson = 0;
int printxml = 0;
int static_config = 0;
int output_depth = FREESASA_OUTPUT_CHAIN;
int no_rel = 0;

// chain groups
int n_chain_groups = 0;
char** chain_groups = NULL;

// selection commands
int n_select = 0;
char** select_cmd = NULL;


void
help(void)
{
    fprintf(stderr, "\nUsage: %s [options] pdb-file(s)\n\n", program_name);
    fprintf(stderr, "GENERAL OPTIONS\n"
            "  -h (--help)           Print this message\n"
            "  -v (--version)        Print version of the program\n");
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
            "  -l (--no-log)         Don't print log message (useful with -r -R and -B)\n"
            "  -w (--no-warnings)    Don't print warnings (will still print warnings due to\n"
            "                        invalid command line options)\n"
            "\n"
            "  -o <file> (--output=<file>)\n"
            "  -e <file> (--error-file=<file>)\n"
            "                        Redirect output and/or errors and warnings to file.\n"
            "\n"
            "  -r  (--foreach-residue-type) or  --residue-type-file=<output-file>\n"
            "  -R  (--foreach-residue)      or  --residue-file=<output-file>\n"
            "                        Print SASA for each residue, either grouped by type or\n"
            "                        sequentially. Use the -file variant to specify an output\n"
            "                        file.\n"
            "\n"
            "  -B  (--print-as-B-values)  or   --B-value-file=<output-file>\n"
            "                        Print PDB file where the temperature factor of each atom\n"
            "                        has been replaced by its SASA, and the occupancy number\n"
            "                        by the atomic radius. Use the -file variant to specify\n"
            "                        an output file.\n"
            "                        This option might give confusing output when used in\n"
        "                        conjuction with the options -C and -g.\n"
        "\n");
    if (USE_JSON) {
        fprintf(stderr,
                "  --json  or  --json-file=<file>\n"
                "                        Print results in JSON format\n"
                "\n");
    }
    if (USE_XML) {
        fprintf(stderr,
                "  --xml   or  --xml-file=<file>\n"
                "                        Print results in XML format\n"
                "\n");
    }
    if (USE_JSON || USE_XML) {
        fprintf(stderr,
                "  --output-depth=<structure|chain|residue|atom>\n"
                "                        Level of detail in JSON and/or XML output [default: chain]\n"
                "\n");
    }
    fprintf(stderr,
            "  --rsa   or  --rsa-file=<file>\n"
            "                        Print relative SASA values in RSA format (same as\n"
            "                        NACCESS). The reference amino acid SASAs are calculated\n"
            "                        using ProtOr radii.\n"
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

void
short_help(void)
{
    fprintf(stderr, "Run '%s -h' for usage instructions.\n",
            program_name);
}

void release_resources()
{
    if (classifier_from_file) freesasa_classifier_free(classifier_from_file);
    if (output_pdb) fclose(output_pdb);
    if (per_residue_type_file) fclose(per_residue_type_file);
    if (per_residue_file) fclose(per_residue_file);
    if (rsa_file) fclose(rsa_file);
    if (json_file) fclose(json_file);
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
    fputc('\n', stderr);
    short_help();
    fputc('\n', stderr);
    fflush(stderr);
    release_resources();
    exit(EXIT_FAILURE);
}

freesasa_structure **
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


void
run_analysis(FILE *input,
             const char *name)
{
    int name_len = strlen(name);
    freesasa_result *result = NULL;
    freesasa_nodearea classes;
    freesasa_structure **structures = NULL;
    freesasa_result_node *tree = freesasa_result_tree_new();
    int n = 0, rel = (no_rel ? FREESASA_OUTPUT_SKIP_REL : 0);

    // read PDB file
    structures = get_structures(input, &n);
    if (n == 0) abort_msg("Invalid input.");
    
    if (printlog) {
        freesasa_write_parameters(output, &parameters);
    }
    
    // perform calculation on each structure and output results
    for (int i = 0; i < n; ++i) {
        char name_i[name_len+10];
        result = freesasa_calc_structure(structures[i], &parameters);
        if (result == NULL)        abort_msg("Can't calculate SASA.");
        classes = freesasa_result_classes(structures[i], result);
        strcpy(name_i,name);
        if (n > 1 && (structure_options & FREESASA_SEPARATE_MODELS))
            sprintf(name_i+strlen(name_i), ":%d", freesasa_structure_model(structures[i]));
        if (printlog) {
            if (n > 1) fprintf(output,"\n\n####################\n");
            freesasa_write_result(output, result, name_i, 
                                  freesasa_structure_chain_labels(structures[i]), &classes);
            freesasa_per_chain(output, result, structures[i]);
        }
        if (per_residue_type) {
            if (n > 1) fprintf(per_residue_type_file, "\n## %s\n", name_i);
            freesasa_per_residue_type(per_residue_type_file, result, structures[i]);
        }
        if (per_residue) {
            if (n > 1) fprintf(per_residue_file, "\n## %s\n", name_i);
            freesasa_per_residue(per_residue_file, result, structures[i]);
        }
        if (printpdb) {
            freesasa_write_pdb(output_pdb, result, structures[i]);
        }
        if (n_select > 0) {
            fprintf(output,"\nSELECTIONS\n");
            for (int c = 0; c < n_select; ++c) {
                double a;
                char sel_name[FREESASA_MAX_SELECTION_NAME+1];
                if (freesasa_select_area(select_cmd[c], sel_name, &a, structures[i], result)
                    == FREESASA_SUCCESS) {
                    fprintf(output, "%s : %10.2f\n", sel_name, a);
                } else {
                    abort_msg("Illegal selection");
                }
            }
        }
        if (printrsa || printjson || printxml) {
            if (freesasa_result_tree_add_result(tree, result, structures[i], name_i)
                != FREESASA_SUCCESS) {
                abort_msg("Error generating result-tree");
            }
        }
        freesasa_result_free(result);
    }
    
    if (printrsa)  freesasa_export_tree(rsa_file,  tree, &parameters, FREESASA_RSA | rel);
    if (printjson) freesasa_export_tree(json_file, tree, &parameters, FREESASA_JSON | output_depth | rel);
    if (printxml)  freesasa_export_tree(xml_file,  tree, &parameters, FREESASA_XML | output_depth | rel);

    freesasa_result_node_free(tree);
    for (int i = 0; i < n; ++i) freesasa_structure_free(structures[i]);
    free(structures);
}

FILE*
fopen_werr(const char* filename,
           const char* mode) 
{
    errno = 0;
    FILE *f = fopen(filename, mode);
    if (f == NULL) {
        abort_msg("could not open file '%s'; %s",
                  filename, strerror(errno));
    }
    return f;
}

void
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

void
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

void
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
    enum {B_FILE, RES_FILE, SEQ_FILE, SELECT, UNKNOWN,
          RSA_FILE, RSA, JSON_FILE, JSON, XML_FILE, XML,
          O_DEPTH, RADII};
    parameters = freesasa_default_parameters;
    memset(opt_set, 0, n_opt);
    program_name = "freesasa";

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
        {"radius-from-occupancy",no_argument,       0, 'O'},
        {"residue-type-file",    required_argument, &option_flag, RES_FILE},
        {"residue-file",         required_argument, &option_flag, SEQ_FILE},
        {"B-value-file",         required_argument, &option_flag, B_FILE},
        {"select",               required_argument, &option_flag, SELECT},
        {"unknown",              required_argument, &option_flag, UNKNOWN},
        {"rsa-file",             required_argument, &option_flag, RSA_FILE},
        {"rsa",                  no_argument,       &option_flag, RSA},
        {"json-file",            required_argument, &option_flag, JSON_FILE},
        {"json",                 no_argument,       &option_flag, JSON},
        {"xml-file",             required_argument, &option_flag, XML_FILE},
        {"xml",                  no_argument,       &option_flag, XML},
        {"output-depth",         required_argument, &option_flag, O_DEPTH},
        {"radii",                required_argument, &option_flag, RADII},
        {0,0,0,0}
    };
    options_string = ":hvlwLSHYOCMmBrRc:n:t:p:g:e:o:";
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
                per_residue_type_file = fopen_werr(optarg, "w");
                break;
            case SEQ_FILE:
                per_residue = 1;
                per_residue_file = fopen_werr(optarg, "w");
                break;
            case B_FILE:
                printpdb = 1;
                output_pdb = fopen_werr(optarg, "w");
                break;
            case SELECT:
                add_select(optarg);
                break;
            case UNKNOWN:
                add_unknown_option(optarg);
                break;
            case RSA:
                printrsa = 1;
                break;
            case RSA_FILE:
                printrsa = 1;
                rsa_file = fopen_werr(optarg, "w");
                break;
            case JSON:
                printjson = 1;
                break;
            case JSON_FILE:
                printjson = 1;
                json_file = fopen_werr(optarg, "w");
                break;
            case XML:
                printxml = 1;
                break;
            case XML_FILE:
                printxml = 1;
                xml_file = fopen_werr(optarg, "w");
                break;
            case O_DEPTH:
                if (strcmp("structure", optarg) == 0)
                    output_depth = FREESASA_OUTPUT_STRUCTURE;
                else if (strcmp("chain", optarg) == 0)
                    output_depth = FREESASA_OUTPUT_CHAIN;
                else if (strcmp("residue", optarg) == 0)
                    output_depth = FREESASA_OUTPUT_RESIDUE;
                else if (strcmp("atom", optarg) == 0)
                    output_depth = 0;
                else abort_msg("Output depth '%s' not allowed, "
                               "can only be 'structure', 'chain', 'residue' or 'atom'",
                               optarg);
                break;
            case RADII:
                static_config = 1;
                if (strcmp("naccess", optarg) == 0) {
                    classifier = &freesasa_naccess_classifier;
                } else if (strcmp("protor", optarg) == 0) {
                    classifier = &freesasa_protor_classifier;
                } else {
                    abort_msg("Config '%s' not allowed, "
                              "can only be 'protor' or 'naccess')", optarg);
                }
                break;
            default:
                abort(); // what does this even mean?
            }
            break;
        case 'h':
            help();
            exit(EXIT_SUCCESS);
        case 'v':
            printf("FreeSASA %s\n", version);
            printf("License: MIT <http://opensource.org/licenses/MIT>\n");
            printf("If you use this program for research, please cite:\n");
            printf("  Simon Mitternacht (2016) FreeSASA: An open source C\n"
                   "  library for solvent accessible surface area calculations.\n"
                   "  F1000Research 5:189.\n");
            exit(EXIT_SUCCESS);
        case 'e': 
            errlog = fopen_werr(optarg, "w");
            freesasa_set_err_out(errlog);
            break;
        case 'o':
            output = fopen_werr(optarg, "w");
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
                abort_msg("error: Resolution needs to be at least 1 (20 recommended minum for S&R, 5 for L&R).");
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
                abort_msg("error: probe radius must be 0 or larger.");
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
            per_residue_type = 1;
            break;
        case 'R':
            per_residue = 1;
            break;
        case 'b':
        case 'B':
            printpdb = 1;
            break;
        case 'g':
            add_chain_groups(optarg);
            break;
        case 't':
#if USE_THREADS
            parameters.n_threads = atoi(optarg);
            if (parameters.n_threads < 1) abort_msg("Number of threads must be 1 or larger.");
#else
            abort_msg("Option '-t' only defined if program compiled with thread support.");
#endif
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
    if (per_residue_type_file == NULL) per_residue_type_file = output;
    if (per_residue_file == NULL) per_residue_file = output;
    if (output_pdb == NULL) output_pdb = output;
    if (rsa_file == NULL) rsa_file = output;
    if (json_file == NULL) json_file = output;
    if (xml_file == NULL) xml_file = output;
    if (alg_set > 1) abort_msg("Multiple algorithms specified.");
    if (opt_set['m'] && opt_set['M']) abort_msg("The options -m and -M can't be combined.");
    if (opt_set['g'] && opt_set['C']) abort_msg("The options -g and -C can't be combined.");
    if (opt_set['c'] && static_config) abort_msg("The options -c and --radii cannot be combined");
    if (opt_set['O'] && static_config) abort_msg("The options -O and --radii cannot be combined");
    if (opt_set['c'] && opt_set['O']) abort_msg("The option -c and -O can't be combined");
    if (printrsa && (opt_set['c'] || opt_set['O'])) {
        freesasa_warn("Will skip REL columns in RSA when custom atomic radii selected.");
    }
    if (printrsa && (opt_set['C'] || opt_set['M']))
        abort_msg("The option --rsa can not be combined with -C or -M. "
                  "The RSA format does not support several results in one file.");
    if (printrsa || printjson || printxml) printlog = 0;
    if (printlog) fprintf(output,"## %s %s ##\n", program_name, version);
    if (argc > optind) {
        for (int i = optind; i < argc; ++i) {
            errno = 0;
            input = fopen_werr(argv[i],"r");
            run_analysis(input, argv[i]);
            fclose(input);
        }
    } else {
        if (!isatty(STDIN_FILENO)) run_analysis(stdin, "stdin");
        else abort_msg("No input.", program_name);
    }

    release_resources();

    return EXIT_SUCCESS;
}
 
