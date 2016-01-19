/*
  Copyright Simon Mitternacht 2013-2016.

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

/**
    This source file contains everything that is in freesasa.h
    interface and did not have a natural home in any of the private
    submodules.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "structure.h"
#include "sasa.h"
#include "pdb.h"
#include "freesasa.h"
#include "classifier.h"
#include "util.h"

#ifdef PACKAGE_VERSION
const char *freesasa_version = PACKAGE_VERSION;
#else
const char *freesasa_version = "";
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

const char *freesasa_alg_names[] = {"Lee & Richards", "Shrake & Rupley"};


// Classifier that classifies each atom according to residue
extern const freesasa_classifier freesasa_residue_classifier;

freesasa_strvp*
freesasa_strvp_new(int n);

freesasa_strvp*
freesasa_result_classify(const freesasa_result *result, 
                         const freesasa_structure *structure,
                         const freesasa_classifier *classifier) 
{
    assert(result);
    assert(structure);

    int n_atoms;
    int n_classes;
    freesasa_strvp *strvp;

    if (classifier == NULL) {
        classifier = &freesasa_default_classifier;
        if (classifier == NULL) return NULL;
    }

    n_atoms = freesasa_structure_n(structure);
    n_classes = classifier->n_classes;
    strvp = freesasa_strvp_new(n_classes+1);
    if (strvp == NULL) {mem_fail(); return NULL;}

    for(int i = 0; i < n_classes; ++i) {
        strvp->string[i] = strdup(classifier->class2str(i,classifier));
        if (strvp->string[i] == NULL) {mem_fail(); return NULL;}
        strvp->value[i] = 0;
    }
    strvp->string[n_classes] = strdup("Unknown");
    if (strvp->string[n_classes] == NULL) {mem_fail(); return NULL;}
    strvp->value[n_classes] = 0;

    for (int i = 0; i < n_atoms; ++i) {
        const char *res_name = freesasa_structure_atom_res_name(structure,i);
        const char *atom_name = freesasa_structure_atom_name(structure,i);
        int c = classifier->sasa_class(res_name,atom_name,classifier);
        if (c == FREESASA_WARN) c = n_classes; // unknown
        strvp->value[c] += result->sasa[i];
    }

    return strvp;
}

void
freesasa_result_free(freesasa_result *r)
{
    if (r) {
        free(r->sasa);
        free(r);
    }
}

static freesasa_result*
freesasa_calc(const coord_t *c, 
              const double *radii,
              const freesasa_parameters *parameters)

{
    assert(c);
    assert(radii);

    freesasa_result *result = malloc(sizeof(freesasa_result));
    int ret;
    const freesasa_parameters *p = parameters;
    
    if (result == NULL) { mem_fail(); return NULL; }    
    if (p == NULL) p = &freesasa_default_parameters;

    result->n_atoms = freesasa_coord_n(c);
    result->sasa = malloc(sizeof(double)*result->n_atoms);
    if(result->sasa == NULL) { mem_fail(); freesasa_result_free(result); return NULL; }

    switch(p->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        ret = freesasa_shrake_rupley(result->sasa, c, radii,
                                     p->probe_radius,
                                     p->shrake_rupley_n_points, 
                                     p->n_threads);
        break;
    case FREESASA_LEE_RICHARDS:
        ret = freesasa_lee_richards(result->sasa, c, radii,
                                    p->probe_radius,
                                    p->lee_richards_n_slices, 
                                    p->n_threads);
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
    if (coord != NULL) result = freesasa_calc(coord,radii,parameters);
    if (coord == NULL || result == NULL) {
       freesasa_result_free(result);
        freesasa_coord_free(coord);
        mem_fail();
        result = NULL;
    }
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
int
freesasa_log(FILE *log,
             freesasa_result *result,
             const char *name,
             const freesasa_parameters *parameters,
             const freesasa_strvp* class_sasa)
{
    assert(log);
    if (freesasa_write_parameters(log,parameters) != FREESASA_FAIL &&
        freesasa_write_result(log,result,name,NULL,class_sasa)
        != FREESASA_FAIL)
        return FREESASA_SUCCESS;
    return fail_msg("");
}

int
freesasa_write_result(FILE *log,
                      freesasa_result *result,
                      const char *name,
                      const char *chains,
                      const freesasa_strvp* class_area)
{
    assert(log);

    fprintf(log,"\nINPUT\n");
    if (name == NULL) fprintf(log,"source  : unknown\n");
    else              fprintf(log,"source  : %s\n",name);
    if (chains != NULL) fprintf(log,"chains  : %s\n",chains);
    fprintf(log,"atoms   : %d\n",result->n_atoms);
    
    fprintf(log,"\nRESULTS (A^2)\n");
    if (class_area == NULL) {
        fprintf(log,"Total : %10.2f\n",result->total);
    } else {
        int m = 6;
        char fmt[21];
        for (int i = 0; i < class_area->n; ++i) {
            int l = strlen(class_area->string[i]);
            m = (l > m) ? l : m;
        }
        sprintf(fmt,"%%-%ds : %%10.2f\n",m);
        fprintf(log,fmt,"Total",result->total);
        for (int i = 0; i < class_area->n; ++i) {
            if (class_area->value[i] > 0)
                fprintf(log,fmt,
                        class_area->string[i],
                        class_area->value[i]);
        }
    }

    fflush(log);
    if (ferror(log)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int freesasa_write_parameters(FILE *log,
                              const freesasa_parameters *parameters)
{
    assert(log);
    const freesasa_parameters *p = parameters;
    if (p == NULL) p = &freesasa_default_parameters;

    fprintf(log,"\nPARAMETERS\n");

    fprintf(log,"algorithm    : %s\n",freesasa_alg_names[p->alg]);
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
                          freesasa_result *result,
                          const freesasa_structure *structure)
{
    assert(output);
    assert(structure);

    const freesasa_classifier *c = &freesasa_residue_classifier;
    freesasa_strvp *residue_area
        = freesasa_result_classify(result,structure,c);

    for (int i = 0; i < c->n_classes; ++i) {
        double sasa = residue_area->value[i];
        if (i < 20 || sasa > 0) {
            fprintf(output,"RES %s : %10.2f\n",
                    residue_area->string[i],sasa);
        }
    }
    freesasa_strvp_free(residue_area);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}


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


int
freesasa_per_residue(FILE *output,
                     freesasa_result *result,
                     const freesasa_structure *structure)
{
    assert(output);
    assert(structure);
    
    const int naa = freesasa_structure_n_residues(structure);
    for (int i = 0; i < naa; ++i) {
        errno = 0;
        int area =
            fprintf(output,"SEQ %s : %7.2f\n",
                    freesasa_structure_residue_descriptor(structure,i),
                    freesasa_single_residue_sasa(result,structure,i));
        if (area < 0)
            return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

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
    for (int i = 0; i < n; ++i) svp->string[i] = NULL;
    return svp;
}

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

int
freesasa_set_verbosity(freesasa_verbosity s) 
{
    if (s == FREESASA_V_NORMAL ||
        s == FREESASA_V_NOWARNINGS ||
        s == FREESASA_V_SILENT) {
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

