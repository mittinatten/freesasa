#if HAVE_CONFIG_H
# include <config.h>
#endif
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#include "pdb.h"
#include "freesasa_internal.h"
#include "classifier.h"

void
freesasa_residue_rel_subarea(freesasa_subarea *rel,
                             const freesasa_subarea *abs,
                             const freesasa_subarea *ref)
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
                 const char *protein_name)
{
#ifdef PACKAGE_VERSION
    fprintf(output, "REM  FreeSASA " PACKAGE_VERSION "\n");
#else
    fprintf(output, "REM  FreeSASA\n");
#endif
    fprintf(output, "REM  Absolute and relative SASAs for %s\n", protein_name);
    fprintf(output, "REM  Reference values calculated using %s radii\n", config_name);
    fprintf(output, "REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar\n");
    fprintf(output, "REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL\n");
}

static inline void 
rsa_print_abs_rel(FILE*output,
                  double abs,
                  double rel)
{
    fprintf(output, "%7.2f", abs);
    if (isfinite(rel)) fprintf(output, "%6.1f", rel);
    else fprintf(output, "   N/A");
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
                  const freesasa_subarea *abs,
                  const freesasa_subarea *rel,
                  const freesasa_structure *structure)
{
    const char *resi_str;
    char chain;

    resi_str = freesasa_structure_residue_number(structure, iaa);
    chain = freesasa_structure_residue_chain(structure, iaa);

    fprintf(output, "RES %s %c%s  ", abs->name, chain, resi_str);
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

int
freesasa_write_rsa(FILE *output,
                   const freesasa_structure_node *tree)
{
    assert(output);
    assert(tree);

    const freesasa_structure *structure = freesasa_structure_node_structure(tree);
    const freesasa_structure_node *residue, *chain = freesasa_structure_node_children(tree);
    const freesasa_subarea *abs, *reference;
    freesasa_subarea rel;
    int res_index, chain_index;

    rsa_print_header(output, freesasa_structure_node_classified_by(tree),
                     freesasa_structure_node_name(tree));

    res_index = chain_index = 0;
    while(chain) {
        residue = freesasa_structure_node_children(chain);
        while (residue) {
            abs = freesasa_structure_node_area(residue);
            reference = freesasa_structure_node_residue_reference(residue);
            if (reference) {
                freesasa_residue_rel_subarea(&rel, abs, reference);
            } else {
                rel = freesasa_subarea_null;
            }
            rsa_print_residue(output, res_index, abs, &rel, structure);
            ++res_index;
            residue = freesasa_structure_node_next(residue);
        }
        chain = freesasa_structure_node_next(chain);
    }

    fprintf(output, "END  Absolute sums over single chains surface\n");

    chain = freesasa_structure_node_children(tree);
    chain_index = 0;
    while(chain) {
        const char *name = freesasa_structure_node_name(chain);
        abs = freesasa_structure_node_area(chain);
        
        fprintf(output,"CHAIN%3d %c %10.1f   %10.1f   %10.1f   %10.1f   %10.1f\n",
                chain_index+1, name[0], abs->total, abs->side_chain,
                abs->main_chain, abs->apolar, abs->polar);

        ++chain_index;
        chain = freesasa_structure_node_next(chain);
    }

    abs = freesasa_structure_node_area(tree);
    fprintf(output, "END  Absolute sums over all chains\n");
    fprintf(output,"TOTAL      %10.1f   %10.1f   %10.1f   %10.1f   %10.1f\n",
            abs->total, abs->side_chain, abs->main_chain, abs->apolar, abs->polar);
    
    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    
    return FREESASA_SUCCESS;
}
