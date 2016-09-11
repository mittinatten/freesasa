/**
    This source file contains everything that is in freesasa.h
    interface and did not have a natural home in any of the private
    submodules.
 */

#if HAVE_CONFIG_H
# include <config.h>
#endif
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <errno.h>

#include "freesasa_internal.h"
#include "pdb.h"
#include "classifier.h"

#ifdef PACKAGE_VERSION
const char *freesasa_version = PACKAGE_VERSION;
#else
const char *freesasa_version = "";
#endif

#ifdef PACKAGE_STRING
const char *freesasa_string = PACKAGE_STRING;
#else
const char *freesasa_string = "FreeSASA";
#endif

// to control error messages (used for debugging and testing)
static freesasa_verbosity verbosity;

// Allows compilation with different defaults
// depending on USE_THREADS. but still exposing the value in a header
// that doesn't depend on USE_THREADS
#ifdef USE_THREADS
#define DEF_NUMBER_THREADS 2
#else
#define DEF_NUMBER_THREADS 1
#endif
const int FREESASA_DEF_NUMBER_THREADS = DEF_NUMBER_THREADS;

const freesasa_parameters freesasa_default_parameters = {
    .alg = FREESASA_DEF_ALGORITHM,
    .probe_radius = FREESASA_DEF_PROBE_RADIUS,
    .shrake_rupley_n_points = FREESASA_DEF_SR_N,
    .lee_richards_n_slices = FREESASA_DEF_LR_N,
    .n_threads = DEF_NUMBER_THREADS,
};

freesasa_strvp*
freesasa_strvp_new(int n);

static freesasa_result *
result_new(int n)
{
    freesasa_result *result = malloc(sizeof(freesasa_result));

    if (result == NULL) {
        mem_fail();
        return NULL;
    }

    result->sasa = malloc(sizeof(double)  * n);

    if (result->sasa == NULL) {
        mem_fail();
        freesasa_result_free(result);
        return NULL;
    }

    result->n_atoms = n;

    return result;
}

void
freesasa_result_free(freesasa_result *r)
{
    if (r) {
        free(r->sasa);
        free(r);
    }
}

freesasa_result*
freesasa_calc(const coord_t *c, 
              const double *radii,
              const freesasa_parameters *parameters)

{
    assert(c);
    assert(radii);

    freesasa_result *result = result_new(freesasa_coord_n(c));
    int ret;

    if (result == NULL) {
        fail_msg("");
        return NULL;
    }

    if (parameters == NULL) parameters = &freesasa_default_parameters;

    switch(parameters->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        ret = freesasa_shrake_rupley(result->sasa, c, radii, parameters);
        break;
    case FREESASA_LEE_RICHARDS:
        ret = freesasa_lee_richards(result->sasa, c, radii, parameters);
        break;
    default:
        assert(0); //should never get here
        break;
    }
    if (ret == FREESASA_FAIL) {
        freesasa_result_free(result);
        return NULL;
    }

    result->total = 0;
    for (int i = 0; i < freesasa_coord_n(c); ++i) {
        result->total += result->sasa[i];
    }
    result->parameters = *parameters;

    return result;
}

freesasa_result*
freesasa_calc_coord(const double *xyz, 
                    const double *radii,
                    int n,
                    const freesasa_parameters *parameters)
{
    assert(xyz);
    assert(radii);
    assert(n > 0);

    coord_t *coord = NULL;
    freesasa_result *result = NULL;

    coord = freesasa_coord_new_linked(xyz,n);
    if (coord != NULL) result = freesasa_calc(coord, radii, parameters);
    if (result == NULL) fail_msg("");
    
    freesasa_coord_free(coord);

    return result;
}

freesasa_result*
freesasa_calc_structure(const freesasa_structure* structure,
                        const freesasa_parameters* parameters)
{
    assert(structure);

    return freesasa_calc(freesasa_structure_xyz(structure),
                         freesasa_structure_radius(structure),
                         parameters);
}

freesasa_result_node *
freesasa_calc_tree(const freesasa_structure *structure,
                   const freesasa_parameters *parameters,
                   const char *name)
{
    assert(structure);

    freesasa_result_node *tree = NULL;
    freesasa_result *result = freesasa_calc(freesasa_structure_xyz(structure),
                                            freesasa_structure_radius(structure),
                                            parameters);

    if (result != NULL) {
        tree = freesasa_result_tree_init(result, structure, name);
    } else {
        fail_msg("");
    }

    if (tree == NULL) {
        fail_msg("");
    }

    freesasa_result_free(result);

    return tree;
}

static int
write_result(FILE *log,
             const freesasa_result_node *result)
{
    assert(log);
    assert(freesasa_result_node_type(result) == FREESASA_NODE_RESULT);

    const char *name = NULL;
    const freesasa_result_node *structure = NULL, *chain = NULL;
    const freesasa_nodearea *area = NULL;

    name = freesasa_result_node_name(result);
    structure = freesasa_result_node_children(result);
    assert(structure);
    area = freesasa_result_node_area(structure);
    assert(area);

    fprintf(log,"\nINPUT\n");
    if (name == NULL) fprintf(log,"source  : unknown\n");
    else              fprintf(log,"source  : %s\n",name);
    fprintf(log,"chains  : %s\n", freesasa_result_node_structure_chain_labels(structure));
    fprintf(log,"atoms   : %d\n", freesasa_result_node_structure_n_atoms(structure));

    fprintf(log,"\nRESULTS (A^2)\n");
    fprintf(log,"Total   : %10.2f\n", area->total);
    fprintf(log,"Apolar  : %10.2f\n", area->apolar);
    fprintf(log,"Polar   : %10.2f\n", area->polar);
    if (area->unknown > 0) {
        fprintf(log,"Unknown : %10.2f\n",area->unknown);
    }

    chain = freesasa_result_node_children(structure);
    while (chain) {
        area = freesasa_result_node_area(chain);
        assert(area);
        fprintf(log, "CHAIN %s : %10.2f\n",freesasa_result_node_name(chain), area->total);
        chain = freesasa_result_node_next(chain);
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

static int
write_results(FILE *log,
              const freesasa_result_node *root)
{
    assert(log);
    assert(freesasa_result_node_type(root) == FREESASA_NODE_ROOT);

    const freesasa_result_node *result = freesasa_result_node_children(root);
    int ret;

    //int several = (freesasa_result_node_next(result) != NULL); // are there more than one result
    while(result) {
        //if (several) fprintf(log, "\n\n####################\n");
        ret = write_result(log, result);
        if (ret != FREESASA_SUCCESS) return ret;
        result = freesasa_result_node_next(result);
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


double
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


int
freesasa_per_residue(FILE *output,
                     const freesasa_result *result,
                     const freesasa_structure *structure)
{
    assert(output);
    assert(structure);

    const int naa = freesasa_structure_n_residues(structure);

    for (int i = 0; i < naa; ++i) {
        int area =
            fprintf(output,"SEQ %c %s %s : %7.2f\n",
                    freesasa_structure_residue_chain(structure, i),
                    freesasa_structure_residue_number(structure, i),
                    freesasa_structure_residue_name(structure, i),
                    freesasa_single_residue_sasa(result, structure, i));
        if (area < 0)
            return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

int
freesasa_export_tree(FILE *file,
                     const freesasa_result_node *root,
                     const freesasa_parameters *parameters,
                     int options)
{
    assert(freesasa_result_node_type(root) == FREESASA_NODE_ROOT);
    if (parameters == NULL) parameters = &freesasa_default_parameters;
    if (options & FREESASA_LOG) {
        return freesasa_write_log(file, root);
    }
    if (options & FREESASA_PDB) {
        return freesasa_write_pdb(file, root);
    }
    if (options & FREESASA_RSA) {
        return freesasa_write_rsa(file, root, parameters, options);
    }
    if (options & FREESASA_JSON) {
        if (USE_JSON) return freesasa_write_json(file, root, parameters, options);
        else return fail_msg("Library was built without support for JSON output.");
    }
    if (options & FREESASA_XML) {
        if (USE_XML) return freesasa_write_xml(file, root, parameters, options);
        else return fail_msg("Library was built without support for XML output.");
    }
    return fail_msg("No valid options given");
}

freesasa_result *
freesasa_result_clone(const freesasa_result *result)
{
    freesasa_result *clone = result_new(result->n_atoms);
    if (clone == NULL) {
        fail_msg("");
        return NULL;
    }

    clone->n_atoms = result->n_atoms;
    clone->total = result->total;
    clone->parameters = result->parameters;
    memcpy(clone->sasa, result->sasa, sizeof(double) * clone->n_atoms);

    return clone;
}

const char*
freesasa_alg_name(freesasa_algorithm alg)
{
    switch(alg) {
    case FREESASA_SHRAKE_RUPLEY:
        return "Shrake & Rupley";
    case FREESASA_LEE_RICHARDS:
        return "Lee & Richards";
    }
    assert(0 && "Illegal algorithm");
}

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

int
freesasa_write_log(FILE *log,
                   const freesasa_result_node *root)
{
    // leave it this way until we have a way of loggin after all
    // calculations are done.
    return write_results(log, root);
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

