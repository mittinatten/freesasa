#if HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <assert.h>
#include "freesasa_internal.h"
#include "classifier.h"

// to control error messages (used for debugging and testing)
static freesasa_verbosity verbosity;

int
freesasa_set_verbosity(freesasa_verbosity s) 
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
    assert(log);
    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    const char *name = NULL;
    freesasa_node *structure = NULL, *chain = NULL;
    const freesasa_nodearea *area = NULL;

    name = freesasa_node_name(result);
    structure = freesasa_node_children(result);
    assert(structure);
    area = freesasa_node_area(structure);
    assert(area);

    fprintf(log,"\nINPUT\n");
    if (name == NULL) fprintf(log,"source  : unknown\n");
    else              fprintf(log,"source  : %s\n",name);
    fprintf(log,"chains  : %s\n", freesasa_node_structure_chain_labels(structure));
    fprintf(log,"atoms   : %d\n", freesasa_node_structure_n_atoms(structure));

    fprintf(log,"\nRESULTS (A^2)\n");
    fprintf(log,"Total   : %10.2f\n", area->total);
    fprintf(log,"Apolar  : %10.2f\n", area->apolar);
    fprintf(log,"Polar   : %10.2f\n", area->polar);
    if (area->unknown > 0) {
        fprintf(log,"Unknown : %10.2f\n",area->unknown);
    }

    chain = freesasa_node_children(structure);
    while (chain) {
        area = freesasa_node_area(chain);
        assert(area);
        fprintf(log, "CHAIN %s : %10.2f\n",freesasa_node_name(chain), area->total);
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
    while (structure) {
        const freesasa_selection **selection = freesasa_node_structure_selections(structure);
        if (selection && *selection) {
            fprintf(log, "\nSELECTIONS\n");
            while(*selection) {
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

int
freesasa_write_parameters(FILE *log,
                          const freesasa_parameters *parameters)
{
    assert(log);
    const freesasa_parameters *p = parameters;
    if (p == NULL) p = &freesasa_default_parameters;

    fprintf(log,"\nPARAMETERS\n");

    fprintf(log,"algorithm    : %s\n",freesasa_alg_name(p->alg));
    fprintf(log,"probe-radius : %.3f\n", p->probe_radius);
    if (USE_THREADS)
        fprintf(log,"threads      : %d\n",p->n_threads);

    switch(p->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        fprintf(log,"testpoints   : %d\n",p->shrake_rupley_n_points);
        break;
    case FREESASA_LEE_RICHARDS:
        fprintf(log,"slices       : %d\n",p->lee_richards_n_slices);
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

int
freesasa_write_res(FILE *log,
                   freesasa_node *root)
{
    assert(log);
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    freesasa_node *result, *structure, *chain, *residue;
    int n_res = freesasa_classify_n_residue_types()+1, i_res;
    double *residue_area = malloc(sizeof(double) * n_res);

    if (residue_area == NULL) return mem_fail();

    result = freesasa_node_children(root);
    while (result) {
        for (int i = 0; i < n_res; ++i) residue_area[i] = 0;
        structure = freesasa_node_children(result);
        while (structure) {
            chain = freesasa_node_children(structure);
            while (chain) {
                residue = freesasa_node_children(chain);
                while(residue) {
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
        for (int i_res = 0; i_res < n_res; ++i_res) {
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

int
freesasa_write_seq(FILE *log,
                   freesasa_node *root)
{
    assert(log);
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);
    freesasa_node *result, *structure, *chain, *residue;

    result = freesasa_node_children(root);
    while (result) {
        structure = freesasa_node_children(result);
        fprintf(log, "# Residues in %s\n", freesasa_node_name(result));
        while (structure) {
            chain = freesasa_node_children(structure);
            while (chain) {
                residue = freesasa_node_children(chain);
                while(residue) {
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
        fprintf(log,"\n");
        result = freesasa_node_next(result);
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int
freesasa_write_log(FILE *log,
                   freesasa_node *root)
{
    assert(log);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    freesasa_node *result = freesasa_node_children(root);
    int several = (freesasa_node_next(result) != NULL); // are there more than one result
    int err = 0;

    while(result) {
        if (several) fprintf(log, "\n\n####################\n");
        if (write_result(log, result) == FREESASA_FAIL) ++err;
        if (write_selections(log, result) == FREESASA_FAIL) ++err;
        result = freesasa_node_next(result);
    }
    if (err) return FREESASA_FAIL;

    return FREESASA_SUCCESS;
}

/** 
 ** Below are a number of output functions that are obsolete with the
 ** node interface / freesasa_export_tree(). They have been marked as
 ** deprecated in the documentation, but could probably just be
 ** deleted in 2.0, since they are not a central part of the API.
 **/
// deprecated
freesasa_strvp*
freesasa_strvp_new(int n);

//deprecated
int
freesasa_per_residue_type(FILE *output,
                          const freesasa_result *result,
                          const freesasa_structure *structure)
{
    assert(output);
    assert(structure);

    int n_res = freesasa_classify_n_residue_types()+1, i_res;
    freesasa_strvp *residue_area = freesasa_strvp_new(n_res);
    if (residue_area == NULL) return FREESASA_FAIL;

    for (int i = 0; i < result->n_atoms; ++i) {
        i_res = freesasa_classify_residue(freesasa_structure_atom_res_name(structure, i));
        residue_area->value[i_res] += result->sasa[i];
    }

    for (int i_res = 0; i_res < n_res; ++i_res) {
        double sasa = residue_area->value[i_res];
        if (i_res < 20 || sasa > 0) {
            fprintf(output,"RES %s : %10.2f\n",
                    freesasa_classify_residue_name(i_res),
                    sasa);
        }
    }
    freesasa_strvp_free(residue_area);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

// deprecated
static double
freesasa_single_residue_sasa(const freesasa_result *r,
                             const freesasa_structure *s,
                             int r_i)
{
    assert(r);
    assert(s);

    int first, last;
    const double *sasa = r->sasa;            
    double a = 0;

    freesasa_structure_residue_atoms(s,r_i,&first,&last);

    for (int j = first; j <= last; ++j) {
        a += sasa[j];
    }

    return a;
}

//deprecated
int
freesasa_per_residue(FILE *output,
                     const freesasa_result *result,
                     const freesasa_structure *structure)
{
    assert(output);
    assert(structure);

    const int naa = freesasa_structure_n_residues(structure);

    for (int i = 0; i < naa; ++i) {
        fprintf(output,"SEQ %c %s %s : %7.2f\n",
                freesasa_structure_residue_chain(structure, i),
                freesasa_structure_residue_number(structure, i),
                freesasa_structure_residue_name(structure, i),
                freesasa_single_residue_sasa(result, structure, i));
    }
    
    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

//deprecated
int
freesasa_per_chain(FILE *output,
                   freesasa_result *result,
                   const freesasa_structure *structure)
{
    const char *chains = freesasa_structure_chain_labels(structure);
    int n_chains = strlen(chains);

    for (int c = 0; c < n_chains; ++c) {
        double area = 0;
        for (int i = 0; i < result->n_atoms; ++i)
            if (freesasa_structure_atom_chain(structure, i) == chains[c])
                area += result->sasa[i];
        fprintf(output,"CHAIN %c : %10.2f\n", chains[c], area);
    }

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

// deprecated
freesasa_strvp*
freesasa_strvp_new(int n)
{
    freesasa_strvp* svp = malloc(sizeof(freesasa_strvp));
    if (svp == NULL) {mem_fail(); return NULL;}
    svp->value = malloc(sizeof(double)*n);
    svp->string = malloc(sizeof(char*)*n);
    if (!svp->value || !svp->string) {
        free(svp->value);  svp->value = NULL;
        free(svp->string); svp->string = NULL;
        freesasa_strvp_free(svp);
        mem_fail();
        return NULL;
    }
    svp->n = n;
    for (int i = 0; i < n; ++i) {
        svp->string[i] = NULL;
        svp->value[i] = 0;
    }
    return svp;
}

//deprecated
freesasa_strvp*
freesasa_result_classify(const freesasa_result *result, 
                         const freesasa_structure *structure,
                         const freesasa_classifier *classifier) 
{
    assert(result);
    assert(structure);

    int n_atoms;
    int n_classes = FREESASA_ATOM_UNKNOWN+1;
    freesasa_strvp *strvp;

    if (classifier == NULL) {
        classifier = &freesasa_default_classifier;
        if (classifier == NULL) return NULL;
    }

    n_atoms = freesasa_structure_n(structure);
    strvp = freesasa_strvp_new(n_classes);
    if (strvp == NULL) {mem_fail(); return NULL;}

    for(int i = 0; i < n_classes; ++i) {
        strvp->string[i] = strdup(freesasa_classifier_class2str(i));
        if (strvp->string[i] == NULL) {mem_fail(); return NULL;}
        strvp->value[i] = 0;
    }

    for (int i = 0; i < n_atoms; ++i) {
        const char *res_name = freesasa_structure_atom_res_name(structure,i);
        const char *atom_name = freesasa_structure_atom_name(structure,i);
        int c = freesasa_classifier_class(classifier, res_name, atom_name);
        strvp->value[c] += result->sasa[i];
    }

    return strvp;
}

//deprecated
void
freesasa_strvp_free(freesasa_strvp *svp)
{
    if (svp) {
        if (svp->value) free(svp->value);
        if (svp->string) {
            for (int i = 0; i < svp->n; ++i) {
                if (svp->string[i]) free(svp->string[i]);
            }
            free(svp->string);
        }
        free(svp);
    }
}
