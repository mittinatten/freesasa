#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>

#include "classifier.h"
#include "freesasa_internal.h"
#include "pdb.h"

void freesasa_residue_rel_nodearea(freesasa_nodearea *rel,
                                   const freesasa_nodearea *abs,
                                   const freesasa_nodearea *ref)
{
    rel->total = 100. * abs->total / ref->total;
    rel->side_chain = 100. * abs->side_chain / ref->side_chain;
    rel->main_chain = 100. * abs->main_chain / ref->main_chain;
    rel->polar = 100. * abs->polar / ref->polar;
    rel->apolar = 100. * abs->apolar / ref->apolar;
    rel->name = abs->name;
}

static void
rsa_print_header(FILE *output,
                 const char *config_name,
                 const char *protein_name,
                 const char *chains,
                 const freesasa_parameters *parameters,
                 int options)
{
    freesasa_algorithm alg = parameters->alg;
#ifdef PACKAGE_VERSION
    fprintf(output, "REM  FreeSASA " PACKAGE_VERSION "\n");
#else
    fprintf(output, "REM  FreeSASA\n");
#endif
    fprintf(output, "REM  Absolute and relative SASAs for %s\n", protein_name);
    if (!(options & FREESASA_OUTPUT_SKIP_REL))
        fprintf(output, "REM  Atomic radii and reference values for relative SASA: %s\n", config_name);
    else
        fprintf(output, "REM  No reference values available to calculate relative SASA\n");
    fprintf(output, "REM  Chains: %s\n", chains);
    fprintf(output, "REM  Algorithm: %s\n", freesasa_alg_name(alg));
    fprintf(output, "REM  Probe-radius: %.2f\n", parameters->probe_radius);
    if (alg == FREESASA_LEE_RICHARDS) {
        fprintf(output, "REM  Slices: %d\n", parameters->lee_richards_n_slices);
    } else if (alg == FREESASA_SHRAKE_RUPLEY) {
        fprintf(output, "REM  Test-points: %d\n", parameters->shrake_rupley_n_points);
    }
    fprintf(output, "REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar\n");
    fprintf(output, "REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL\n");
}

static inline void
rsa_print_abs_rel(FILE *output,
                  double abs,
                  double rel)
{
    fprintf(output, "%7.2f", abs);
    if (isfinite(rel))
        fprintf(output, "%6.1f", rel);
    else
        fprintf(output, "   N/A");
}

static inline void
rsa_print_abs_only(FILE *output,
                   double abs)
{
    fprintf(output, "%7.2f", abs);
    fprintf(output, "   N/A");
}

static int
rsa_print_residue(FILE *output,
                  int iaa,
                  const freesasa_nodearea *abs,
                  const freesasa_nodearea *rel,
                  freesasa_node *residue)
{
    const char *resi_str;
    char chain;

    resi_str = freesasa_node_residue_number(residue);
    chain = freesasa_node_name(freesasa_node_parent(residue))[0];

    fprintf(output, "RES %s %c%s ", abs->name, chain, resi_str);
    if (rel->name != NULL) {
        rsa_print_abs_rel(output, abs->total, rel->total);
        rsa_print_abs_rel(output, abs->side_chain, rel->side_chain);
        rsa_print_abs_rel(output, abs->main_chain, rel->main_chain);
        rsa_print_abs_rel(output, abs->apolar, rel->apolar);
        rsa_print_abs_rel(output, abs->polar, rel->polar);
    } else {
        rsa_print_abs_only(output, abs->total);
        rsa_print_abs_only(output, abs->side_chain);
        rsa_print_abs_only(output, abs->main_chain);
        rsa_print_abs_only(output, abs->apolar);
        rsa_print_abs_only(output, abs->polar);
    }
    fprintf(output, "\n");
    return FREESASA_SUCCESS;
}

int freesasa_write_rsa(FILE *output,
                       freesasa_node *tree,
                       int options)
{
    freesasa_node *residue, *chain, *structure_node, *result_node;
    const freesasa_nodearea *abs, *reference;
    freesasa_nodearea rel;
    int res_index, chain_index;
    const freesasa_parameters *parameters;

    assert(output);
    assert(tree);

    result_node = freesasa_node_children(tree);
    parameters = freesasa_node_result_parameters(result_node);
    structure_node = freesasa_node_children(result_node);
    chain = freesasa_node_children(structure_node);

    rsa_print_header(output, freesasa_node_classified_by(result_node),
                     freesasa_node_name(result_node), freesasa_node_name(structure_node), parameters, options);

    res_index = chain_index = 0;
    while (chain) {
        residue = freesasa_node_children(chain);
        while (residue) {
            abs = freesasa_node_area(residue);
            reference = freesasa_node_residue_reference(residue);
            if (reference && !(options & FREESASA_OUTPUT_SKIP_REL)) {
                freesasa_residue_rel_nodearea(&rel, abs, reference);
            } else {
                rel = freesasa_nodearea_null;
            }
            rsa_print_residue(output, res_index, abs, &rel, residue);
            ++res_index;
            residue = freesasa_node_next(residue);
        }
        chain = freesasa_node_next(chain);
    }

    fprintf(output, "END  Absolute sums over single chains surface\n");

    chain = freesasa_node_children(structure_node);
    chain_index = 0;
    while (chain) {
        const char *name = freesasa_node_name(chain);
        abs = freesasa_node_area(chain);

        fprintf(output, "CHAIN%3d %c %10.1f   %10.1f   %10.1f   %10.1f   %10.1f\n",
                chain_index + 1, name[0], abs->total, abs->side_chain,
                abs->main_chain, abs->apolar, abs->polar);

        ++chain_index;
        chain = freesasa_node_next(chain);
    }

    abs = freesasa_node_area(structure_node);
    fprintf(output, "END  Absolute sums over all chains\n");
    fprintf(output, "TOTAL      %10.1f   %10.1f   %10.1f   %10.1f   %10.1f\n",
            abs->total, abs->side_chain, abs->main_chain, abs->apolar, abs->polar);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}
