#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <errno.h>
#include <stdlib.h>

#include "classifier.h"
#include "freesasa_internal.h"

/** to control error messages (used for debugging and testing) */
static freesasa_verbosity verbosity;

int freesasa_set_verbosity(freesasa_verbosity s)
{
    if (s == FREESASA_V_NORMAL ||
        s == FREESASA_V_NOWARNINGS ||
        s == FREESASA_V_SILENT ||
        s == FREESASA_V_DEBUG) {
        verbosity = s;
        return FREESASA_SUCCESS;
    }
    return FREESASA_WARN;
}

freesasa_verbosity
freesasa_get_verbosity(void)
{
    return verbosity;
}

static int
write_result(FILE *log,
             freesasa_node *result)
{
    const char *name = NULL;
    freesasa_node *structure = NULL, *chain = NULL;
    const freesasa_nodearea *area = NULL;

    assert(log);
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    name = freesasa_node_name(result);
    structure = freesasa_node_children(result);
    assert(structure);
    area = freesasa_node_area(structure);
    assert(area);

    fprintf(log, "\nINPUT\n");
    if (name == NULL)
        fprintf(log, "source  : unknown\n");
    else
        fprintf(log, "source  : %s\n", name);
    fprintf(log, "chains  : %s\n", freesasa_node_structure_chain_labels(structure));
    fprintf(log, "model   : %d\n", freesasa_node_structure_model(structure));
    fprintf(log, "atoms   : %d\n", freesasa_node_structure_n_atoms(structure));

    fprintf(log, "\nRESULTS (A^2)\n");
    fprintf(log, "Total   : %10.2f\n", area->total);
    fprintf(log, "Apolar  : %10.2f\n", area->apolar);
    fprintf(log, "Polar   : %10.2f\n", area->polar);
    if (area->unknown > 0) {
        fprintf(log, "Unknown : %10.2f\n", area->unknown);
    }

    chain = freesasa_node_children(structure);
    while (chain) {
        area = freesasa_node_area(chain);
        assert(area);
        fprintf(log, "CHAIN %s : %10.2f\n", freesasa_node_name(chain), area->total);
        chain = freesasa_node_next(chain);
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

static int
write_selections(FILE *log,
                 freesasa_node *result)
{
    freesasa_node *structure = freesasa_node_children(result);
    const freesasa_selection **selection;

    while (structure) {
        selection = freesasa_node_structure_selections(structure);
        if (selection && *selection) {
            fprintf(log, "\nSELECTIONS\n");
            while (*selection) {
                fprintf(log, "%s : %10.2f\n",
                        freesasa_selection_name(*selection),
                        freesasa_selection_area(*selection));
                ++selection;
            }
        }
        structure = freesasa_node_next(structure);
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

static int
write_parameters(FILE *log,
                 const freesasa_parameters *parameters)
{
    const freesasa_parameters *p = parameters;

    assert(log);

    if (p == NULL) p = &freesasa_default_parameters;

    fprintf(log, "\nPARAMETERS\n");

    fprintf(log, "algorithm    : %s\n", freesasa_alg_name(p->alg));
    fprintf(log, "probe-radius : %.3f\n", p->probe_radius);
#if USE_THREADS
    fprintf(log, "threads      : %d\n", p->n_threads);
#endif

    switch (p->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        fprintf(log, "testpoints   : %d\n", p->shrake_rupley_n_points);
        break;
    case FREESASA_LEE_RICHARDS:
        fprintf(log, "slices       : %d\n", p->lee_richards_n_slices);
        break;
    default:
        assert(0);
        break;
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int freesasa_write_res(FILE *log,
                       freesasa_node *root)
{
    freesasa_node *result, *structure, *chain, *residue;
    int n_res = freesasa_classify_n_residue_types() + 1, i_res, i;
    double *residue_area = malloc(sizeof(double) * n_res);

    assert(log);
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    if (residue_area == NULL) return mem_fail();

    result = freesasa_node_children(root);
    while (result) {
        for (i = 0; i < n_res; ++i)
            residue_area[i] = 0;
        structure = freesasa_node_children(result);
        while (structure) {
            chain = freesasa_node_children(structure);
            while (chain) {
                residue = freesasa_node_children(chain);
                while (residue) {
                    assert(freesasa_node_type(residue) == FREESASA_NODE_RESIDUE);
                    i_res = freesasa_classify_residue(freesasa_node_name(residue));
                    residue_area[i_res] += freesasa_node_area(residue)->total;
                    residue = freesasa_node_next(residue);
                }
                chain = freesasa_node_next(chain);
            }
            structure = freesasa_node_next(structure);
        }

        fprintf(log, "# Residue types in %s\n", freesasa_node_name(result));
        for (i_res = 0; i_res < n_res; ++i_res) {
            double sasa = residue_area[i_res];
            if (i_res < 20 || sasa > 0) {
                fprintf(log, "RES %s : %10.2f\n",
                        freesasa_classify_residue_name(i_res),
                        sasa);
            }
        }
        fprintf(log, "\n");

        result = freesasa_node_next(result);
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int freesasa_write_seq(FILE *log,
                       freesasa_node *root)
{
    freesasa_node *result, *structure, *chain, *residue;

    assert(log);
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    result = freesasa_node_children(root);

    while (result) {
        structure = freesasa_node_children(result);
        fprintf(log, "# Residues in %s\n", freesasa_node_name(result));
        while (structure) {
            chain = freesasa_node_children(structure);
            while (chain) {
                residue = freesasa_node_children(chain);
                while (residue) {
                    assert(freesasa_node_type(residue) == FREESASA_NODE_RESIDUE);
                    fprintf(log, "SEQ %s %s %s : %7.2f\n",
                            freesasa_node_name(chain),
                            freesasa_node_residue_number(residue),
                            freesasa_node_name(residue),
                            freesasa_node_area(residue)->total);
                    residue = freesasa_node_next(residue);
                }
                chain = freesasa_node_next(chain);
            }
            structure = freesasa_node_next(structure);
        }
        fprintf(log, "\n");
        result = freesasa_node_next(result);
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int freesasa_write_log(FILE *log,
                       freesasa_node *root)
{
    freesasa_node *result = freesasa_node_children(root);
    int several = (freesasa_node_next(result) != NULL); /* are there more than one result */
    int err = 0;

    assert(log);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    if (write_parameters(log, freesasa_node_result_parameters(result)) == FREESASA_FAIL)
        ++err;

    while (result) {
        if (several) fprintf(log, "\n\n####################\n");
        if (write_result(log, result) == FREESASA_FAIL) ++err;
        if (write_selections(log, result) == FREESASA_FAIL) ++err;
        result = freesasa_node_next(result);
    }
    if (err) return FREESASA_FAIL;

    return FREESASA_SUCCESS;
}
