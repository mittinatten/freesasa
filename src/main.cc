#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <iostream>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "cif.hh"
#include "freesasa.h"

#if STDC_HEADERS
extern int getopt(int, char *const *, const char *);
extern int optind, optopt;
extern char *optarg;
#endif

static const char *program_name = "freesasa";

#if USE_XML
#define XML_STRING "|xml"
#else
#define XML_STRING ""
#endif

#if USE_JSON
#define JSON_STRING "|json"
#else
#define JSON_STRING ""
#endif

#define FORMAT_STRING "log|res|seq|pdb|rsa|cif" XML_STRING JSON_STRING

enum { B_FILE,
       SELECT,
       UNKNOWN,
       RSA,
       RADII,
       DEPRECATED,
       CIF };

static int option_flag;

static struct option long_options[] = {
    {"lee-richards", no_argument, 0, 'L'},
    {"shrake-rupley", no_argument, 0, 'S'},
    {"probe-radius", required_argument, 0, 'p'},
    {"resolution", required_argument, 0, 'n'},
    {"help", no_argument, 0, 'h'},
    {"version", no_argument, 0, 'v'},
    {"no-warnings", no_argument, 0, 'w'},
    {"n-threads", required_argument, 0, 't'},
    {"config-file", required_argument, 0, 'c'},
    {"radius-from-occupancy", no_argument, 0, 'O'},
    {"hetatm", no_argument, 0, 'H'},
    {"hydrogen", no_argument, 0, 'Y'},
    {"separate-chains", no_argument, 0, 'C'},
    {"separate-models", no_argument, 0, 'M'},
    {"join-models", no_argument, 0, 'm'},
    {"chain-groups", required_argument, 0, 'g'},
    {"error-file", required_argument, 0, 'e'},
    {"output", required_argument, 0, 'o'},
    {"format", required_argument, 0, 'f'},
    {"depth", required_argument, 0, 'd'},
    {"cif", no_argument, &option_flag, CIF},
    {"select", required_argument, &option_flag, SELECT},
    {"unknown", required_argument, &option_flag, UNKNOWN},
    {"rsa", no_argument, &option_flag, RSA},
    {"radii", required_argument, &option_flag, RADII},
    {"deprecated", no_argument, &option_flag, DEPRECATED},
    /* Deprecated options */
    {"foreach-residue-type", no_argument, 0, 'r'},
    {"foreach-residue", no_argument, 0, 'R'},
    {"print-as-B-values", no_argument, 0, 'B'},
    {"no-log", no_argument, 0, 'l'},
    {0, 0, 0, 0}};

#define NOARG_OPTIONS "hvwLSHYOCMm"
#define NOARG_DEPRECATED "BrRl"
#define ARG_OPTIONS "c:n:t:p:g:e:o:f:"
const char *options_string = ":" NOARG_OPTIONS NOARG_DEPRECATED ARG_OPTIONS;

/* State of app (most settings are stored here) */
struct cli_state {
    freesasa_parameters parameters;
    freesasa_classifier *classifier_from_file;
    const freesasa_classifier *classifier;
    int structure_options;
    int static_classifier;
    int cif;
    int no_rel;
    /* chain groups */
    int n_chain_groups;
    char **chain_groups;
    /* selection commands */
    int n_select;
    char **select_cmd;
    /* output settings */
    int output_format, output_depth;
    /* Files */
    char *output_filename;
    FILE *input, *output, *errlog;
};

struct analysis_results {
    freesasa_node *tree;
    freesasa_structure **structures;
    int n_structures;
};

/* Defaults */
static void
init_state(struct cli_state *state)
{
    state->parameters = freesasa_default_parameters;
    state->classifier_from_file = NULL;
    state->classifier = NULL;
    state->structure_options = 0;
    state->static_classifier = 0;
    state->no_rel = 0;
    state->n_chain_groups = 0;
    state->chain_groups = NULL;
    state->n_select = 0;
    state->select_cmd = 0;
    state->output_format = 0;
    state->output_depth = FREESASA_OUTPUT_CHAIN;
    state->output_filename = NULL;
    state->output = NULL;
    state->errlog = NULL;
    state->cif = 0;
}

static void
release_state(struct cli_state *state)
{
    int i;

    if (state->classifier_from_file)
        freesasa_classifier_free(state->classifier_from_file);
    if (state->chain_groups) {
        for (i = 0; i < state->n_chain_groups; ++i) {
            free(state->chain_groups[i]);
        }
    }
    if (state->select_cmd) {
        for (i = 0; i < state->n_select; ++i) {
            free(state->select_cmd[i]);
        }
    }
    if (state->errlog) fclose(state->errlog);
    if (state->output) fclose(state->output);
    free(state->output_filename);
}

static void
addresses(FILE *out)
{
    fprintf(out, "\n" REPORTBUG "\nHome page: " HOMEPAGE "\n");
}

static void
help(void)
{
    printf("\nUsage: %s [options] pdb-file ...", program_name);
    printf("\n       %s [options] < pdb-file", program_name);
    printf("\n       %s (--help | --version | --deprecated)\n", program_name);
    printf("\n"
           "Options:\n"
           "  --shrake-rupley | --lee-richards\n"
           "  --probe-radius=<NUMBER>\n"
           "  --resolution=<INTEGER> -n-threads=<INTEGER>\n"
           "  --radius-from-occupancy | --config-file=<FILE> | --radii=<protor|naccess>\n"
           "  --hetatm --hydrogen\n"
           "  --unknown=<guess|skip|halt>\n"
           "  --cif\n"
           "  --separate-models | --join-models\n"
           "  --separate-chains | --chain-groups=<LIST> ...\n"
           "  --select=<STRING> ...\n"
           "  --output=<FILE> --error-file=<FILE> --no-warnings\n"
           "  --format=<" FORMAT_STRING "> ... \n"
           "  --depth=<structure|chain|residue|atom>\n");
    printf("\nPlease refer to the man pages or online documentation for more information.\n");
    addresses(stdout);
    printf("\n");
}

static void
deprecated(void)
{
    fprintf(stderr,
            "These options will disappear in later versions of FreeSASA.\n"
            "Use --format instead\n\n"
            "  --rsa                         Equivalent to --format=rsa\n"
            "  -B  --print-as-B-values       Equivalent to --format=pdb\n"
            "  -r  --foreach-residue-type    Equivalent to --format=res\n"
            "  -R  --foreach-residue         Equivalent to --format=seq.\n"
            "  -l  --no-log                  Log suppressed if other format selected.\n"
            "                                Option has no effect."
            "\n");
}

static void
short_help(void)
{
    fprintf(stderr, "Call '%s -h' or 'man %s' for usage instructions\n",
            program_name, program_name);
    addresses(stderr);
}

static void
version(void)
{
    printf("%s\n", PACKAGE_STRING);
    printf("License: MIT <http://opensource.org/licenses/MIT>\n");
    printf("If you use this program for research, please cite:\n");
    printf("  Simon Mitternacht (2016) FreeSASA: An open source C\n"
           "  library for solvent accessible surface area calculations.\n"
           "  F1000Research 5:189.\n");
    addresses(stdout);
}

/** error handling **/

static void
err_msg(const char *prefix,
        const char *format,
        ...)
{
    va_list arg;
    va_start(arg, format);
    fprintf(stderr, "%s: %s: ", program_name, prefix);
    vfprintf(stderr, format, arg);
    fputc('\n', stderr);
    va_end(arg);
    fflush(stderr);
}

static void
exit_with_help(void)
{
    fputc('\n', stderr);
    short_help();
    fputc('\n', stderr);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

#define warn(...) err_msg("warning", __VA_ARGS__)
#define error(...) err_msg("error", __VA_ARGS__)
#define abort_msg(...)      \
    do {                    \
        error(__VA_ARGS__); \
        exit_with_help();   \
    } while (0)

static std::vector<freesasa_structure *>
get_structures(std::FILE *input,
               const struct cli_state *state)
{
    int n = 0, n2 = 0;
    std::vector<freesasa_structure *> structures;
    freesasa_structure *tmp;

    if ((state->structure_options & FREESASA_SEPARATE_CHAINS) ||
        (state->structure_options & FREESASA_SEPARATE_MODELS)) {
        if (state->cif) {
            structures = freesasa_cif_structure_array(input, &n, state->classifier, state->structure_options);
        } else {
            // TODO this hack needed since PDB implementation is in C
            freesasa_structure **db_ptr_structs = freesasa_structure_array(input, &n, state->classifier, state->structure_options);
            structures.reserve(n);
            for (int i = 0; i < n; ++i) {
                structures.push_back(std::move(db_ptr_structs[i]));
            }
        }
    } else {
        n = 1;
        if (state->cif) {
            structures.emplace_back(freesasa_structure_from_cif(input, state->classifier, state->structure_options));
        } else {
            structures.emplace_back(freesasa_structure_from_pdb(input, state->classifier, state->structure_options));
        }
        if (structures[0] == NULL) {
            abort_msg("invalid input");
        }
    }

    /* get chain-groups (if requested) */
    if (state->n_chain_groups > 0) {
        n2 = n;
        for (int i = 0; i < state->n_chain_groups; ++i) {
            for (int j = 0; j < n; ++j) {
                // TODO make this function pdb and cif compatible.
                tmp = freesasa_structure_get_chains(structures[j], state->chain_groups[i],
                                                    state->classifier, state->structure_options);
                if (tmp != NULL) {
                    ++n2;
                    structures.reserve(n2);
                    structures.push_back(tmp);
                } else {
                    abort_msg("at least one of chain(s) '%s' not found", state->chain_groups[i]);
                }
            }
        }
        n = n2;
    }

    if (n == 0) {
        abort_msg("invalid input");
    }

    return structures;
}

static freesasa_node *
run_analysis(std::vector<freesasa_structure *> structures,
             const char *name,
             const struct cli_state *state)
{
    int name_len = strlen(name);
    freesasa_node *tree = freesasa_tree_new(), *tmp_tree, *structure_node;
    const freesasa_result *result;
    freesasa_selection *sel;
    char *name_i = (char *)malloc(name_len + 10);

    if (tree == NULL) abort_msg("failed to initialize result-tree");
    if (name_i == NULL) abort_msg("memory failure");

    /* perform calculation on each structure */
    for (auto structure : structures) {
        strcpy(name_i, name);
        if (structures.size() > 1 && (state->structure_options & FREESASA_SEPARATE_MODELS))
            snprintf(name_i + strlen(name_i), 9, ":%d", freesasa_structure_model(structure));

        tmp_tree = freesasa_calc_tree(structure, &state->parameters, name_i);
        if (tmp_tree == NULL) abort_msg("can't calculate SASA");

        structure_node =
            freesasa_node_children(freesasa_node_children(tmp_tree));
        result = freesasa_node_structure_result(structure_node);

        /* Calculate selections for each structure */
        if (state->n_select > 0) {
            for (int c = 0; c < state->n_select; ++c) {
                sel = freesasa_selection_new(state->select_cmd[c], structure, result);
                if (sel != NULL) {
                    freesasa_node_structure_add_selection(structure_node, sel);
                } else {
                    abort_msg("illegal selection");
                }
                freesasa_selection_free(sel);
            }
        }

        if (freesasa_tree_join(tree, &tmp_tree) != FREESASA_SUCCESS) {
            abort_msg("failed joining result-trees");
        }
    }

    return tree;
}

static FILE *
fopen_werr(const char *filename,
           const char *mode)
{
    FILE *f = fopen(filename, mode);

    if (f == NULL) {
        abort_msg("could not open file '%s'; %s",
                  filename, strerror(errno));
    }
    return f;
}

static void
state_add_chain_groups(const char *cmd, struct cli_state *state)
{
    int err = 0;
    char *str;
    const char *token;
    char a;

    /* check that string is valid */
    for (size_t i = 0; i < strlen(cmd); ++i) {
        a = cmd[i];
        if (a != '+' &&
            !(a >= 'a' && a <= 'z') && !(a >= 'A' && a <= 'Z') &&
            !(a >= '0' && a <= '9')) {
            error("character '%c' not valid chain ID in --chain-groups, "
                  "valid characters are [A-z0-9] and '+' as separator",
                  a);
            ++err;
        }
    }

    /* extract chain groups */
    if (err == 0) {
        str = strdup(cmd);
        token = strtok(str, "+");
        while (token) {
            ++state->n_chain_groups;
            state->chain_groups = (char **)realloc(state->chain_groups, sizeof(char *) * state->n_chain_groups);
            if (state->chain_groups == NULL) {
                abort_msg("out of memory");
            }
            state->chain_groups[state->n_chain_groups - 1] = strdup(token);
            token = strtok(0, "+");
        }
        free(str);
    } else {
        abort_msg("aborting");
    }
}

static void
state_add_select(const char *cmd, struct cli_state *state)
{
    ++state->n_select;
    state->select_cmd = (char **)realloc(state->select_cmd, sizeof(char *) * state->n_select);
    if (state->select_cmd == NULL) {
        abort_msg("out of memory");
    }
    state->select_cmd[state->n_select - 1] = strdup(cmd);
    if (state->select_cmd[state->n_select - 1] == NULL) {
        abort_msg("out of memory");
    }
}

static void
state_add_unknown_option(const char *optarg, struct cli_state *state)
{
    if (strcmp(optarg, "skip") == 0) {
        state->structure_options |= FREESASA_SKIP_UNKNOWN;
        return;
    }
    if (strcmp(optarg, "halt") == 0) {
        state->structure_options |= FREESASA_HALT_AT_UNKNOWN;
        return;
    }
    if (strcmp(optarg, "guess") == 0) {
        return; /* default */
    }
    abort_msg("unknown alternative to option --unknown: '%s'", optarg);
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
        if (USE_JSON) {
            return FREESASA_JSON;
        } else {
            abort_msg("program was built without JSON support");
        }
    }
    if (strcmp(optarg, "xml") == 0) {
        if (USE_XML) {
            return FREESASA_XML;
        } else {
            abort_msg("program was built without XML support");
        }
    }
    if (strcmp(optarg, "pdb") == 0) {
        return FREESASA_PDB;
    }
    if (strcmp(optarg, "cif") == 0) {
        return FREESASA_CIF;
    }
    abort_msg("unknown output format: '%s'", optarg);
    return FREESASA_FAIL; /* to avoid compiler warnings */
}

static int
parse_output_depth(const char *optarg)
{
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
    abort_msg("output depth '%s' not allowed, "
              "can only be 'structure', 'chain', 'residue' or 'atom'",
              optarg);
    return FREESASA_FAIL; /* to avoid compiler warnings */
}

static void
state_set_static_classifier(const char *optarg, struct cli_state *state)
{
    if (strcmp("naccess", optarg) == 0) {
        state->classifier = &freesasa_naccess_classifier;
    } else if (strcmp("protor", optarg) == 0) {
        state->classifier = &freesasa_protor_classifier;
    } else {
        abort_msg("config '%s' not allowed, "
                  "can only be 'protor' or 'naccess')",
                  optarg);
    }
    state->static_classifier = 1;
}

/* Parse command line arguments and transform state
   accordingly. Parameter state assumed to be initialized to default. */
static int
parse_arg(int argc, char **argv, struct cli_state *state)
{
    int alg_set = 0;
    int opt;
    int n_opt = 'z' + 1;
    char opt_set['z' + 1];
    int option_index = 0;
    FILE *cf;

    memset(opt_set, 0, n_opt);

    while ((opt = getopt_long(argc, argv, options_string,
                              long_options, &option_index)) != -1) {
        opt_set[opt] = 1;
        /* Assume arguments starting with dash are actually missing arguments */
        if (optarg != NULL && optarg[0] == '-') {
            if (option_index > 0)
                abort_msg("missing argument? Value '%s' cannot be argument to '--%s'.\n",
                          program_name, optarg, long_options[option_index].name);
            else
                abort_msg("missing argument? Value '%s' cannot be argument to '-%c'.\n",
                          optarg, opt);
        }
        switch (opt) {
        case 0:
            switch (long_options[option_index].val) {
            case SELECT:
                state_add_select(optarg, state);
                break;
            case UNKNOWN:
                state_add_unknown_option(optarg, state);
                break;
            case RSA:
                state->output_format = FREESASA_RSA;
                break;
            case RADII:
                state_set_static_classifier(optarg, state);
                break;
            case DEPRECATED:
                deprecated();
                exit(EXIT_SUCCESS);
            case CIF:
                state->cif = 1;
                break;
            default:
                abort(); /* what does this even mean? */
            }
            break;
        case 'h':
            help();
            exit(EXIT_SUCCESS);
        case 'v':
            version();
            exit(EXIT_SUCCESS);
        case 'e':
            state->errlog = fopen_werr(optarg, "w");
            freesasa_set_err_out(state->errlog);
            break;
        case 'o':
            if (state->output_filename != NULL) {
                abort_msg("option --output can only be set once");
            }
            state->output_filename = strdup(optarg);
            break;
        case 'f':
            state->output_format |= parse_output_format(optarg);
            break;
        case 'd':
            state->output_depth = parse_output_depth(optarg);
            break;
        case 'w':
            freesasa_set_verbosity(FREESASA_V_NOWARNINGS);
            break;
        case 'c': {
            cf = fopen_werr(optarg, "r");
            state->classifier = state->classifier_from_file = freesasa_classifier_from_file(cf);
            if (state->classifier_from_file == NULL) abort_msg("can't read file '%s'", optarg);
            state->no_rel = 1;
            break;
        }
        case 'n':
            state->parameters.shrake_rupley_n_points = atoi(optarg);
            state->parameters.lee_richards_n_slices = atoi(optarg);
            if (state->parameters.shrake_rupley_n_points <= 0)
                abort_msg("resolution needs to be at least 1 (20 recommended minum for S&R, 5 for L&R)");
            break;
        case 'S':
            state->parameters.alg = FREESASA_SHRAKE_RUPLEY;
            ++alg_set;
            break;
        case 'L':
            state->parameters.alg = FREESASA_LEE_RICHARDS;
            ++alg_set;
            break;
        case 'p':
            state->parameters.probe_radius = atof(optarg);
            if (state->parameters.probe_radius <= 0)
                abort_msg("probe radius must be 0 or larger");
            break;
        case 'H':
            state->structure_options |= FREESASA_INCLUDE_HETATM;
            break;
        case 'Y':
            state->structure_options |= FREESASA_INCLUDE_HYDROGEN;
            break;
        case 'O':
            state->structure_options |= FREESASA_RADIUS_FROM_OCCUPANCY;
            state->no_rel = 1;
            break;
        case 'M':
            state->structure_options |= FREESASA_SEPARATE_MODELS;
            break;
        case 'm':
            state->structure_options |= FREESASA_JOIN_MODELS;
            break;
        case 'C':
            state->structure_options |= FREESASA_SEPARATE_CHAINS;
            break;
        case 'g':
            state_add_chain_groups(optarg, state);
            break;
        case 't':
            if (USE_THREADS) {
                state->parameters.n_threads = atoi(optarg);
                if (state->parameters.n_threads < 1) abort_msg("number of threads must be 1 or larger");
            } else {
                abort_msg("option '-t' only defined if program compiled with thread support");
            }
            break;
        /* Deprecated options */
        case 'r':
            warn("option '-r' deprecated, use '-f res' or '--format=res' instead");
            state->output_format |= FREESASA_RES;
            break;
        case 'R':
            warn("option '-R' deprecated, use '-f seq' or '--format=seq' instead");
            state->output_format |= FREESASA_SEQ;
            break;
        case 'B':
            warn("option '-B' deprecated, use '-f pdb' or '--format=pdb' instead");
            state->output_format |= FREESASA_PDB;
            break;
        case 'l':
            warn("option '-l' deprecated, has no effect.");
            break;
        /* Errors */
        case ':':
            abort_msg("option '-%c' missing argument");
            break;
        case '?':
        default:
            if (optopt == 0) {
                abort_msg("unknown option '%s'", argv[optind - 1]);
            } else {
                abort_msg("unknown option '-%c'", optopt);
            }
            break;
        }
    }

    if (state->output_filename) {
        state->output = fopen_werr(state->output_filename, "w");
    } else {
        state->output = stdout;
    }

    if (alg_set > 1) abort_msg("multiple algorithms specified");
    if (state->output_format == 0) state->output_format = FREESASA_LOG;
    if (opt_set['m'] && opt_set['M']) abort_msg("the options -m and -M can't be combined");
    if (opt_set['g'] && opt_set['C']) abort_msg("the options -g and -C can't be combined");
    if (opt_set['c'] && state->static_classifier) abort_msg("the options -c and --radii cannot be combined");
    if (opt_set['O'] && state->static_classifier) abort_msg("the options -O and --radii cannot be combined");
    if (opt_set['c'] && opt_set['O']) abort_msg("the options -c and -O can't be combined");
    if (state->output_format == FREESASA_RSA && (opt_set['c'] || opt_set['O'])) {
        warn("will skip REL columns in RSA when custom atomic radii selected");
    }
    if (state->output_format == FREESASA_RSA && (opt_set['C'] || opt_set['M']))
        abort_msg("the RSA format can not be used with the options -C or -M, "
                  "it does not support several results in one file");

    if (state->output_format & FREESASA_LOG) {
        fprintf(state->output, "## %s ##\n", PACKAGE_STRING);
    }

    if (state->output_format == FREESASA_CIF && state->cif != 1) abort_msg("CIF output can not be generated from .pdb input");
    if (state->output_format == FREESASA_PDB && state->cif == 1) abort_msg("PDB output can not be generated from .cif input.");
    if (state->output_format == FREESASA_PDB &&
        state->structure_options & FREESASA_SEPARATE_CHAINS &&
        state->structure_options & FREESASA_SEPARATE_MODELS)
        abort_msg("Cannot output a cif/pdb file with both --separate-chains and --separate-models set. Pick one.");

    return optind;
}

int main(int argc,
         char **argv)
{
    struct cli_state state;
    FILE *input = NULL;
    int optind = 0, ret;

    freesasa_node *tree = freesasa_tree_new(), *tmp;
    if (tree == NULL) abort_msg("error initializing calculation");

    init_state(&state);

    optind = parse_arg(argc, argv, &state);
    std::vector<freesasa_structure *> structures;

    if (argc > optind) {
        for (int i = optind; i < argc; ++i) {
            input = fopen_werr(argv[i], "r");
            structures = get_structures(input, &state);
            tmp = run_analysis(structures, argv[i], &state);
            freesasa_tree_join(tree, &tmp);
            fclose(input);
        }
    } else {
        if (!isatty(STDIN_FILENO)) {
            structures = get_structures(stdin, &state);
            tmp = run_analysis(structures, "stdin", &state);
            freesasa_tree_join(tree, &tmp);
        } else {
            abort_msg("no input", program_name);
        }
    }

    if (state.output_format & FREESASA_CIF) {
        ret = freesasa_export_tree_to_cif(state.output_filename, tree);
    } else {
        ret = freesasa_tree_export(state.output, tree, state.output_format | state.output_depth | (state.no_rel ? FREESASA_OUTPUT_SKIP_REL : 0));
    }
    freesasa_node_free(tree);

    for (auto structure : structures) {
        freesasa_structure_free(structure);
    }

    release_state(&state);

    if (ret == FREESASA_FAIL) {
        freesasa_fail("Error when writing results");
    }

    return ret != FREESASA_FAIL ? EXIT_SUCCESS : EXIT_FAILURE;
}
