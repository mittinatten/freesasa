/*
  Copyright Simon Mitternacht 2013.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <sys/time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "sasa.h"
#include "pdb.h"
#include "sasalib.h"
#include "srp.h"
#include "classify.h"

#define NBUF 100

const char *sasalib_name = "sasalib";

typedef struct {
    int *class;
    int *residue;
} class_t;

typedef struct {
    double *sasa;
    double total;
    double *class;
    double *residue;
} result_t;

struct sasalib_ {
    sasalib_algorithm alg;
    double *r;
    sasalib_coord_t *coord;
    class_t *class;
    int n_atoms;
    double probe_radius;
    int n_sr;
    double d_lr;
    int n_threads;
    double elapsed_time;
    int owns_r;
    int calculated;
    result_t *result;
    char proteinname[SASALIB_NAME_LIMIT];
}; 

const sasalib_t sasalib_def_param = {
    .alg = SASALIB_SHRAKE_RUPLEY,
    .r = NULL,
    .coord = NULL,
    .class = NULL,
    .result = NULL,
    .n_atoms = 0,
    .probe_radius = SASALIB_DEF_PROBE_RADIUS,
    .n_sr = SASALIB_DEF_SR_N,
    .d_lr = SASALIB_DEF_LR_D,
    .n_threads = 1,
    .elapsed_time = 0,
    .owns_r = 0,
    .calculated = 0,
    .proteinname = "undef"
};

const char *sasalib_alg_names[] = {"Lee & Richards", "Shrake & Rupley"};

static class_t* sasalib_classify_structure(const sasalib_structure_t *p)
{
    assert(p != NULL);
    assert(sasalib_structure_n(p) >= 0 && 
	   "Trying to classify atoms of illegal structure.");
    
    size_t n = sasalib_structure_n(p);
    
    class_t* c = (class_t*)malloc(sizeof(class_t));
    c->class = (int*)malloc(sizeof(int)*n);
    c->residue = (int*)malloc(sizeof(int)*n);

    for (int i = 0; i < n; ++i) {
	const char *res_name = sasalib_structure_atom_res_name(p,i);
	const char *atom_name = sasalib_structure_atom_name(p,i);
	c->class[i] = sasalib_classify_class(res_name,atom_name);
	c->residue[i] = sasalib_classify_residue(res_name);
    }

    return c;
}

static void sasalib_class_free(class_t *c) 
{
    if (! c) return;
    free(c->class);
    free(c->residue);
    free(c);
}

static result_t* sasalib_result_new(size_t n_atoms)
{
    result_t *r = (result_t*)malloc(sizeof(result_t));

    int nc = sasalib_classify_nclasses();
    int nr = sasalib_classify_nresiduetypes();
    
    r->sasa = (double*)malloc(sizeof(double)*n_atoms);
    r->class = (double*)malloc(sizeof(double)*nc);
    r->residue = (double*)malloc(sizeof(double)*nr);

    for (int i = 0; i < nc; ++i) r->class[i] = 0;
    for (int i = 0; i < nr; ++i) r->residue[i] = 0;
    
    return r;
}

static void sasalib_result_free(result_t *r)
{
    if (! r) return;
    free(r->sasa);
    free(r->class);
    free(r->residue);
    free(r);
}

static void sasalib_get_class_result(sasalib_t *s, sasalib_structure_t *p)
{
    if (s->class) sasalib_class_free(s->class);
    s->class = sasalib_classify_structure(p);
    
    for (size_t i = 0; i < sasalib_structure_n(p); ++i) {
        s->result->class[s->class->class[i]] += s->result->sasa[i];
	s->result->residue[s->class->residue[i]] += s->result->sasa[i];
    }
}

int sasalib_calc(sasalib_t *s, const sasalib_coord_t *c, const double *r)
{
    struct timeval t1, t2;
    int res = 0;

    gettimeofday(&t1,NULL);
    
    if (s->result) sasalib_result_free(s->result);
    result_t *result;
    s->result = result = sasalib_result_new(s->n_atoms);
    
    switch(s->alg) {
    case SASALIB_SHRAKE_RUPLEY:
        res = sasalib_shrake_rupley(result->sasa, c, r, 
				    s->probe_radius,
				    s->n_sr, s->n_threads);
        break;
    case SASALIB_LEE_RICHARDS:
        res = sasalib_lee_richards(result->sasa, c, r, 
				   s->probe_radius,
				   s->d_lr, s->n_threads);
	break;
    default:
        fprintf(stderr,"%s: error: no SASA algorithm specified.\n",
		sasalib_name);
        return SASALIB_FAIL;
    }
    result->total = 0;
    for (int i = 0; i < sasalib_coord_n(c); ++i) {
	result->total += result->sasa[i];
    }

    gettimeofday(&t2,NULL);
    
    s->elapsed_time = (t2.tv_sec-t1.tv_sec); 
    s->elapsed_time += (t2.tv_usec-t1.tv_usec) / 1e6; // s
    
    s->calculated = 1;

    return res;
}

sasalib_t* sasalib_init()
{
    sasalib_t *s = (sasalib_t*) malloc(sizeof(sasalib_t));
    *s = sasalib_def_param;
    return s;
}

void sasalib_copy_param(sasalib_t *target, const sasalib_t *source)
{
    target->alg = source->alg;
    target->n_sr = source->n_sr;
    target->d_lr = source->d_lr;
    target->n_threads = source->n_threads;
    target->probe_radius = source->probe_radius;
}

void sasalib_free(sasalib_t *s)
{
    if (! s) return;
    if (s->owns_r) {
        free(s->r);
    }
    if (s->class) sasalib_class_free(s->class);
    if (s->result) sasalib_result_free(s->result);
    if (s->coord) sasalib_coord_free(s->coord);
    free(s);
}

int sasalib_calc_coord(sasalib_t *s, const double *coord, 
                       const double *r, size_t n)
{
    s->n_atoms = n;
    sasalib_coord_t *c = sasalib_coord_new_linked(coord,n);
    int res = sasalib_calc(s,c,r);
    sasalib_coord_free(c);
    return res;
}

int sasalib_calc_pdb(sasalib_t *s, FILE *pdb_file)
{
    sasalib_structure_t *p = sasalib_structure_init_from_pdb(pdb_file);
    if (!p) {
	fprintf(stderr,"%s: error: failed to read pdb-file.\n",
		sasalib_name);
	sasalib_structure_free(p);
	return SASALIB_FAIL;
    }
    if (!(s->n_atoms = sasalib_structure_n(p))) {
	fprintf(stderr,"%s: error: pdb-file was empty.\n",
		sasalib_name);
	sasalib_structure_free(p);
	return SASALIB_FAIL;
    }
    s->n_atoms = sasalib_structure_n(p);

    //calc OONS radii
    if (s->r) free(s->r);
    s->r = (double*) malloc(sizeof(double)*sasalib_structure_n(p));
    s->owns_r = 1;
    sasalib_structure_r_def(s->r,p);

    int res = sasalib_calc(s,sasalib_structure_xyz(p),s->r);
    if (!res) sasalib_get_class_result(s,p);
    sasalib_structure_free(p);
    return res;
}

int sasalib_link_coord(sasalib_t *s, const double *coord,
                       double *r, size_t n) 
{
    s->coord = sasalib_coord_new_linked(coord,n);

    if (s->r) free(s->r);
    s->r = r;
    s->owns_r = 0;

    return SASALIB_SUCCESS;
}

int refresh(sasalib_t *s)
{
    if (! s->coord || ! s->r ) {
        fprintf(stderr,"%s: error: trying to refresh unitialized "
                "sasalib_t-object.\n",sasalib_name);
        return SASALIB_FAIL;
    }
    sasalib_calc(s,s->coord,s->r);
    return SASALIB_SUCCESS;
}

int sasalib_set_algorithm(sasalib_t *s, sasalib_algorithm alg)
{
    if (alg == SASALIB_SHRAKE_RUPLEY || alg == SASALIB_LEE_RICHARDS) { 
        s->alg = alg;
        return SASALIB_SUCCESS;
    }
    fprintf(stderr,"%s: warning: undefined algorithm selected, "
            "proceeding with default.\n",sasalib_name);
    return SASALIB_WARN;
}

sasalib_algorithm sasalib_get_algorithm(const sasalib_t *s)
{
    return s->alg;
}

const char* sasalib_algorithm_name(const sasalib_t *s)
{
    assert(s->alg == SASALIB_SHRAKE_RUPLEY || s->alg == SASALIB_LEE_RICHARDS);
    return sasalib_alg_names[s->alg];
}

int sasalib_set_probe_radius(sasalib_t *s,double r) 
{
    if (r < 0 || !isfinite(r)) {
	fprintf(stderr,"%s: error: Probe radius r = %f not allowed.",
		sasalib_name, r);
	return SASALIB_WARN;
    }
    s->probe_radius = r;
    return SASALIB_SUCCESS;
}

double sasalib_get_probe_radius(const sasalib_t *s)
{
    return s->probe_radius;
}

int sasalib_set_sr_points(sasalib_t *s, int n) {
    if (sasalib_srp_n_is_valid(n)) { 
        s->n_sr = n;
        return SASALIB_SUCCESS;
    }
    fprintf(stderr,"%s: error: Number of test-points must be"
	    " one of the following:  ", sasalib_name);
    sasalib_srp_print_n_opt(stderr);
    fprintf(stderr,"%s: error: Proceeding with default value (%d)\n",
            sasalib_name,SASALIB_DEF_SR_N);
    s->n_sr = SASALIB_DEF_SR_N;
    return SASALIB_FAIL;
}

int sasalib_get_sr_points(const sasalib_t* s)
{
    if (s->alg == SASALIB_SHRAKE_RUPLEY) return s->n_sr;
    return SASALIB_WARN;
}

int sasalib_set_lr_delta(sasalib_t *s, double d)
{
    if (d > 0 && d <= 5) {
        s->d_lr = d;
        return SASALIB_SUCCESS;
    }
    fprintf(stderr,"%s: error: slice width has to lie between 0 and 5 Å\n",
	    sasalib_name);
    fprintf(stderr,"%s: error: Proceeding with default value (%f)\n",
            sasalib_name, SASALIB_DEF_LR_D);
    s->d_lr = SASALIB_DEF_LR_D;
    return SASALIB_WARN;
}

double sasalib_get_lr_delta(const sasalib_t *s)
{
    if (s->alg == SASALIB_LEE_RICHARDS) return s->d_lr;
    return -1.0;
}

#ifdef PTHREADS
int sasalib_set_nthreads(sasalib_t *s,int n)
{
    if ( n <= 0) {
        fprintf(stderr,"%s: error: Number of threads has to be positive.\n",
		sasalib_name);
	fprintf(stderr,"%s: error: proceeding with default value (1).\n",
		sasalib_name);
        return SASALIB_WARN;
    }
    s->n_threads = n;
    return SASALIB_SUCCESS;
}

int sasalib_get_nthreads(const sasalib_t *s)
{
    return s->n_threads;
}
#endif

size_t sasalib_n_atoms(const sasalib_t *s)
{
    return s->n_atoms;
}

double sasalib_area_total(const sasalib_t *s)
{
    if (! s->calculated) {
        fprintf(stderr,"%s: error: SASA calculation has not been performed, "
                "no total SASA value available.\n",sasalib_name);
        return -1.0;
    }
    return s->result->total;
}
double sasalib_area_class(const sasalib_t* s, sasalib_class c)
{    
    if (! s->calculated) {
        fprintf(stderr,"%s: error: SASA calculation has not been performed, "
                "no SASA value available.\n",sasalib_name);
	return -1.0;
    }
    assert(c >= SASALIB_POLAR && c <= SASALIB_CLASS_UNKNOWN &&
	   "Invalid arguments to sasalib_area_class(2)");
    return s->result->class[c];
}

const double* sasalib_area_atom_array(const sasalib_t *s)
{
    if (s->result->sasa == NULL) {
        fprintf(stderr,"%s: error: SASA calculation has not been performed, "
                "no atomic SASA values are available.\n",sasalib_name);
        return NULL;
    }
    return s->result->sasa;
}
double sasalib_radius_atom(const sasalib_t *s, int i)
{
    if (s->r == NULL) {
	fprintf(stderr,"%s: error: No atomic radii have been assigned.\n",
		sasalib_name);
	return -1.0;
    }
    if (i < 0 || i > s->n_atoms) {
	fprintf(stderr,"%s: error: Atom index %d invalid.\n",
		sasalib_name,i);
	return -1.0;
    }
    return s->r[i];
}
const double* sasalib_radius_atom_array(const sasalib_t *s)
{
    if (s->r == NULL) {
	fprintf(stderr,"%s: error: No atomic radii have been assigned.\n",
		sasalib_name);
	return NULL;
    }
    return s->r;
}
void sasalib_set_proteinname(sasalib_t *s,const char *name)
{
    int n;
    if ((n = strlen(name)) > SASALIB_NAME_LIMIT) {
        strcpy(s->proteinname,"...");
        sprintf(s->proteinname+3,"%.*s",SASALIB_NAME_LIMIT-3,
                name+n+3-SASALIB_NAME_LIMIT);
    } else {
        strcpy(s->proteinname,name);
    }
}

const char* sasalib_get_proteinname(const sasalib_t *s) 
{
    return s->proteinname;
}

int sasalib_log(FILE *log, const sasalib_t *s)
{
    if (! s->calculated) {
        fprintf(stderr, "%s: warning: sasalib_log(2) called, but no calculation"
                "has been performed.\n", sasalib_name);
        fprintf(log, "error: sasalib_log(2) called, but no calculation"
                "has been performed.\n");
        return SASALIB_WARN;
    }
    fprintf(log,"# Using van der Waals radii and atom classes defined \n"
            "# by Ooi et al (PNAS 1987, 84:3086-3090) and a probe radius\n"
            "# of %f Å.\n\n", s->probe_radius);
    fprintf(log,"name: %s\n",s->proteinname);
    fprintf(log,"algorithm: %s\n",sasalib_alg_names[s->alg]);
#ifdef PTHREADS
    fprintf(log,"n_thread: %d\n",s->n_threads);
#endif
    
    switch(s->alg) {
    case SASALIB_SHRAKE_RUPLEY:
        fprintf(log,"n_testpoint: %d\n",s->n_sr);
        break;
    case SASALIB_LEE_RICHARDS:
        fprintf(log,"d_slice: %f Å\n",s->d_lr);
        break;
    default:
        break;
    }
    fprintf(log,"time_elapsed: %f s\n",s->elapsed_time);
    fprintf(log,"n_atoms: %d\n", s->n_atoms);
    return SASALIB_SUCCESS;
}

int sasalib_per_residue(FILE *output, const sasalib_t *s)
{
    if (! s->calculated) {
        fprintf(stderr,"%s: warning: sasalib_per_residue(2) called, "
                "but no calculation has been performed\n",sasalib_name);
        return SASALIB_WARN;
    }
    for (int i = 0; i < sasalib_classify_nresiduetypes(); ++i) {
        double sasa = s->result->residue[i];
        if (i < 20 || sasa > 0) {
            fprintf(output,"%s %8.2f\n",sasalib_classify_residue2str(i),sasa);
        }
    }
    return SASALIB_SUCCESS;
}

