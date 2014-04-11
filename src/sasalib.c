/*
  Copyright Simon Mitternacht 2013-2014.

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
#include <stdarg.h>
#include <assert.h>
#include <math.h>
#if HAVE_CONFIG_H
# include <config.h>
#endif

#include "structure.h"
#include "sasa.h"
#include "pdb.h"
#include "sasalib.h"
#include "srp.h"
#include "classify.h"

#define NBUF 100
#define DEF_NTHREADS 1

#ifdef PACKAGE_NAME
const char *sasalib_name = PACKAGE_NAME;
#else
const char *sasalib_name = "sasalib";
#endif

// to control error messages (used for debugging and testing)
static int sasalib_verbosity;

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
    double *r;
    sasalib_coord_t *coord;
    class_t *class;
    int n_atoms;
    sasalib_structure_t *structure;
    
    //parameters
    sasalib_algorithm alg;
    double probe_radius;
    int n_sr;
    double d_lr;
    int n_threads;
    char proteinname[SASALIB_NAME_LIMIT];
    
    //some internal flags
    int owns_r;
    int calculated;
    
    //results
    double elapsed_time;
    result_t *result;
}; 

const sasalib_t sasalib_def_param = {
    .r = NULL,
    .coord = NULL,
    .class = NULL,
    .n_atoms = 0,
    .structure = NULL,

    .alg = SASALIB_SHRAKE_RUPLEY,
    .probe_radius = SASALIB_DEF_PROBE_RADIUS,
    .n_sr = SASALIB_DEF_SR_N,
    .d_lr = SASALIB_DEF_LR_D,
    .n_threads = DEF_NTHREADS,
    .proteinname = "undef",

    .owns_r = 0,
    .calculated = 0,

    .elapsed_time = 0,
    .result = NULL
};

const char *sasalib_alg_names[] = {"Lee & Richards", "Shrake & Rupley"};

// The following functions lets test-suite silence error messages
enum {SASALIB_V_NORMAL = 0, SASALIB_V_SILENT=1};

int sasalib_set_verbosity(int s) {
    if (s == SASALIB_V_NORMAL || s == SASALIB_V_SILENT) {
	sasalib_verbosity = s;
	return SASALIB_SUCCESS;
    }
    return SASALIB_WARN;
}
int sasalib_get_verbosity() {
    return sasalib_verbosity;
}

static void sasalib_err_impl(int err, const char *format, va_list arg)
{
    fprintf(stderr, "%s: ", sasalib_name);
    switch (err) {
    case SASALIB_FAIL: fputs("error: ", stderr); break; 
    case SASALIB_WARN: fputs("warning: ", stderr); break;
    default: break;
    }
    vfprintf(stderr, format, arg);
    va_end(arg);
    fputc('\n', stderr);
}

int sasalib_fail(const char *format,...) 
{
    if (sasalib_verbosity == SASALIB_V_SILENT) return SASALIB_FAIL;
    va_list arg;
    va_start(arg, format);
    sasalib_err_impl(SASALIB_FAIL,format,arg);
    va_end(arg);
    return SASALIB_FAIL;
}

int sasalib_warn(const char *format,...) 
{
    if (sasalib_verbosity == SASALIB_V_SILENT) return SASALIB_WARN;
    va_list arg;
    va_start(arg, format);
    sasalib_err_impl(SASALIB_WARN,format,arg);
    va_end(arg);
    return SASALIB_WARN;
}

static class_t* sasalib_classify_structure(const sasalib_structure_t *p)
{
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

static int sasalib_calc(sasalib_t *s, 
			const sasalib_coord_t *c, const double *r)
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
        return sasalib_fail("no SASA algorithm specified.");
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
    target->calculated = 0;
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
    if (s->structure) sasalib_structure_free(s->structure);
    free(s);
}

int sasalib_calc_coord(sasalib_t *s, const double *coord, 
                       const double *r, size_t n)
{
    s->calculated = 0;
    s->n_atoms = n;
    s->r = r;

    sasalib_coord_t *c = sasalib_coord_new_linked(coord,n);
    int res = sasalib_calc(s,c,r);
    sasalib_coord_free(c);
    return res;
}

int sasalib_calc_pdb(sasalib_t *s, FILE *pdb_file)
{
    s->calculated = 0;
    sasalib_structure_t *p = sasalib_structure_init_from_pdb(pdb_file);
    if (!p) {
        return sasalib_fail("Failure reading PDB-file.");
    }
    s->n_atoms = sasalib_structure_n(p);
    
    if (s->structure) sasalib_structure_free(s->structure);
    s->structure = p;

    //calc OONS radii
    if (s->r) free(s->r);
    s->r = (double*) malloc(sizeof(double)*sasalib_structure_n(p));
    s->owns_r = 1;
    sasalib_structure_r_def(s->r,p);

    int res = sasalib_calc(s,sasalib_structure_xyz(p),s->r);
    if (!res) sasalib_get_class_result(s,p);

    return res;
}
double sasalib_radius(const char* residue_name, const char* atom_name)
{
    return sasalib_classify_radius(residue_name,atom_name);
}
int sasalib_link_coord(sasalib_t *s, const double *coord,
                       double *r, size_t n) 
{
    s->calculated = 0;
    if (s->coord) sasalib_coord_free(s->coord);
    s->coord = sasalib_coord_new_linked(coord,n);

    if (s->r && s->owns_r) free(s->r);
    s->r = r;
    s->owns_r = 0;

    if (s->result) {
	sasalib_result_free(s->result);
	s->result = NULL;
    }

    s->n_atoms = n;

    return SASALIB_SUCCESS;
}

int sasalib_refresh(sasalib_t *s)
{
    s->calculated = 0;
    if (! s->coord ) 
        return sasalib_fail("sasalib_refresh(1) called, but no coordinates available.");
    if (! s->r )
        return sasalib_fail("sasalib_refresh(1) called, but no atomic radii specified.");
    return sasalib_calc(s,s->coord,s->r);
}

int sasalib_set_algorithm(sasalib_t *s, sasalib_algorithm alg)
{
    s->calculated = 0;
    if (alg == SASALIB_SHRAKE_RUPLEY || alg == SASALIB_LEE_RICHARDS) { 
        s->alg = alg;
        return SASALIB_SUCCESS;
    }
    return sasalib_warn("undefined algorithm selected, proceeding with previously selected algorithm (or default).");
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
    s->calculated = 0;
    if (r < 0 || !isfinite(r)) {
        return sasalib_warn("Probe radius r = %f not allowed.",r);
    }
    s->probe_radius = r;
    return SASALIB_SUCCESS;
}

double sasalib_get_probe_radius(const sasalib_t *s)
{
    return s->probe_radius;
}

int sasalib_set_sr_points(sasalib_t *s, int n) {
    s->calculated = 0;
    if (sasalib_srp_n_is_valid(n)) { 
        s->n_sr = n;
        return SASALIB_SUCCESS;
    }
    sasalib_warn("number of test-points must be one of the following:  ");
    if (sasalib_verbosity == SASALIB_V_NORMAL) {
	sasalib_srp_print_n_opt(stderr);
    }
    sasalib_warn("proceeding with default number of test-points (%d)",SASALIB_DEF_SR_N);
    s->n_sr = SASALIB_DEF_SR_N;
    return SASALIB_WARN;
}

int sasalib_get_sr_points(const sasalib_t* s)
{
    if (s->alg == SASALIB_SHRAKE_RUPLEY) return s->n_sr;
    return SASALIB_WARN;
}

int sasalib_set_lr_delta(sasalib_t *s, double d)
{
    s->calculated = 0;
    if (d > 0 && isfinite(d)) {
        s->d_lr = d;
        if (d <= 2) return SASALIB_SUCCESS;
        if (d > 2) {
            sasalib_warn("for regular SASA calculations, " 
                         "slice width should be less than 2 A.");
            sasalib_warn("proceeding with selected value (%f), "
                         "but results might be inaccurate.", d);
            return SASALIB_WARN;
	}
    }
    sasalib_warn("slice width %f invalid.",d);
    sasalib_warn("proceeding with default value (%f).",SASALIB_DEF_LR_D);
    s->d_lr = SASALIB_DEF_LR_D;
    return SASALIB_WARN;
}

double sasalib_get_lr_delta(const sasalib_t *s)
{
    if (s->alg == SASALIB_LEE_RICHARDS) return s->d_lr;
    return -1.0;
}

#if HAVE_LIBPTHREAD
int sasalib_set_nthreads(sasalib_t *s,int n)
{
    if ( n <= 0) {
	s->n_threads = DEF_NTHREADS;
        return sasalib_warn("Number of threads has to be positive. "
			   "Proceeding with default number (%d).", DEF_NTHREADS);
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
        sasalib_fail("SASA calculation has not been performed, "
                     "no total SASA value available.");
        return -1.0;
    }
    return s->result->total;
}
double sasalib_area_class(const sasalib_t* s, sasalib_class c)
{    
    if (! s->calculated) {
        sasalib_fail("SASA calculation has not been performed, "
                     "no total SASA value available.");
	return -1.0;
    }
    assert(c >= SASALIB_POLAR && c <= SASALIB_CLASS_UNKNOWN &&
           "Invalid arguments to sasalib_area_class(2)");
    return s->result->class[c];
}
double sasalib_area_atom(const sasalib_t *s, int i) 
{
    if ( !s->calculated ) {
	sasalib_fail("SASA calculation has not been performed, "
                     "no SASA values available.");
	return -1.0;
    }
    if (i < 0 || i > s->n_atoms-1) {
        sasalib_fail("Atom index %d invalid.",i);
        return -1.0;
    }
    return s->result->sasa[i];
}

const double* sasalib_area_atom_array(const sasalib_t *s)
{
    if ( !s->calculated ) {
        sasalib_fail("SASA calculation has not been performed, "
                     "no SASA values available.");
        return NULL;
    }
    return s->result->sasa;
}
double sasalib_radius_atom(const sasalib_t *s, int i)
{
    if (s->r == NULL) {
        sasalib_fail("No atomic radii have been assigned.");
        return -1.0;
    }
    if (i < 0 || i > s->n_atoms-1) {
        sasalib_fail("Atom index %d invalid.",i);
        return -1.0;
    }
    return s->r[i];
}
const double* sasalib_radius_atom_array(const sasalib_t *s)
{
    if (s->r == NULL) {
        sasalib_fail("No atomic radii have been assigned.");
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
        const char *msg = "sasalib_log(2) called, but no calculation "
            "has been performed.";
        sasalib_warn(msg);
	if (sasalib_verbosity == SASALIB_V_NORMAL) {
	    fprintf(log,"%s\n",msg);
	}
        return SASALIB_WARN;
    }
    fprintf(log,"name: %s\n",s->proteinname);
    fprintf(log,"algorithm: %s\n",sasalib_alg_names[s->alg]);
    fprintf(log,"probe-radius: %f A\n",s->probe_radius);
#if HAVE_LIBPTHREAD
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
    if (! output) {
	return sasalib_fail("sasalib_per_residue(2) output file is "
			    "a NULL pointer.");
    }
    if (! s->calculated) {
        return sasalib_fail("sasalib_per_residue(2) called, "
                            "but no cqalculation has been performed.");
    }
    for (int i = 0; i < sasalib_classify_nresiduetypes(); ++i) {
        double sasa = s->result->residue[i];
        if (i < 20 || sasa > 0) {
            fprintf(output,"%s %8.2f\n",sasalib_classify_residue2str(i),sasa);
        }
    }
    return SASALIB_SUCCESS;
}

double sasalib_area_residue(const sasalib_t *s, const char *res_name)
{
    if (! s->calculated) {
        sasalib_warn("sasalib_area_residue(2) called, "
                     "but no calculation has been performed.");
        return -1.0;
    }
    int res = sasalib_classify_residue(res_name);
    return s->result->residue[res];
}

int sasalib_write_pdb(FILE *output, const sasalib_t *s)
{
    if (!s->calculated) {
        return sasalib_fail("Cannot output results before "
                            "calculation has been performed.");
    }
    if (!s->structure) {
        return sasalib_fail("Cannot write B-factors: "
                            "no input structure.");
    }
    if (sasalib_structure_n(s->structure) != s->n_atoms) {
        return sasalib_fail("Cannot write B-factors: "
                            "number of atoms mismatch.");
    }
    if (!s->result->sasa) {
        return sasalib_fail("Cannot write B-factors: "
                            "no SASA values available.");
    }
    return sasalib_structure_write_pdb_bfactors(output,	s->structure,
                                                s->result->sasa);
}
