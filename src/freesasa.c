/*
  Copyright Simon Mitternacht 2013-2014.

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
#include "freesasa.h"
#include "srp.h"
#include "classify.h"

#define NBUF 100
#define DEF_NTHREADS 1

#ifdef PACKAGE_NAME
const char *freesasa_name = PACKAGE_NAME;
#else
const char *freesasa_name = "freesasa";
#endif

// to control error messages (used for debugging and testing)
static int freesasa_verbosity;

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

struct freesasa_ {
    double *r;
    freesasa_coord_t *coord;
    class_t *class;
    int n_atoms;
    freesasa_structure_t *structure;
    
    //parameters
    freesasa_algorithm alg;
    double probe_radius;
    int n_sr;
    double d_lr;
    int n_threads;
    char proteinname[FREESASA_NAME_LIMIT];
    
    //some internal flags
    int owns_r;
    int calculated;
    
    //results
    double elapsed_time;
    result_t *result;
}; 

const freesasa_t freesasa_def_param = {
    .r = NULL,
    .coord = NULL,
    .class = NULL,
    .n_atoms = 0,
    .structure = NULL,

    .alg = FREESASA_SHRAKE_RUPLEY,
    .probe_radius = FREESASA_DEF_PROBE_RADIUS,
    .n_sr = FREESASA_DEF_SR_N,
    .d_lr = FREESASA_DEF_LR_D,
    .n_threads = DEF_NTHREADS,
    .proteinname = "undef",

    .owns_r = 0,
    .calculated = 0,

    .elapsed_time = 0,
    .result = NULL
};

const char *freesasa_alg_names[] = {"Lee & Richards", "Shrake & Rupley"};

// The following functions lets test-suite silence error messages
enum {FREESASA_V_NORMAL = 0, FREESASA_V_SILENT=1};

int freesasa_set_verbosity(int s) {
    if (s == FREESASA_V_NORMAL || s == FREESASA_V_SILENT) {
	freesasa_verbosity = s;
	return FREESASA_SUCCESS;
    }
    return FREESASA_WARN;
}
int freesasa_get_verbosity() {
    return freesasa_verbosity;
}

static void freesasa_err_impl(int err, const char *format, va_list arg)
{
    fprintf(stderr, "%s: ", freesasa_name);
    switch (err) {
    case FREESASA_FAIL: fputs("error: ", stderr); break; 
    case FREESASA_WARN: fputs("warning: ", stderr); break;
    default: break;
    }
    vfprintf(stderr, format, arg);
    va_end(arg);
    fputc('\n', stderr);
}

int freesasa_fail(const char *format,...) 
{
    if (freesasa_verbosity == FREESASA_V_SILENT) return FREESASA_FAIL;
    va_list arg;
    va_start(arg, format);
    freesasa_err_impl(FREESASA_FAIL,format,arg);
    va_end(arg);
    return FREESASA_FAIL;
}

int freesasa_warn(const char *format,...) 
{
    if (freesasa_verbosity == FREESASA_V_SILENT) return FREESASA_WARN;
    va_list arg;
    va_start(arg, format);
    freesasa_err_impl(FREESASA_WARN,format,arg);
    va_end(arg);
    return FREESASA_WARN;
}

static class_t* freesasa_classify_structure(const freesasa_structure_t *p)
{
    assert(freesasa_structure_n(p) >= 0 && 
	   "Trying to classify atoms of illegal structure.");
    
    size_t n = freesasa_structure_n(p);
    
    class_t* c = (class_t*)malloc(sizeof(class_t));
    c->class = (int*)malloc(sizeof(int)*n);
    c->residue = (int*)malloc(sizeof(int)*n);

    for (int i = 0; i < n; ++i) {
	const char *res_name = freesasa_structure_atom_res_name(p,i);
	const char *atom_name = freesasa_structure_atom_name(p,i);
	c->class[i] = freesasa_classify_class(res_name,atom_name);
	c->residue[i] = freesasa_classify_residue(res_name);
    }

    return c;
}

static void freesasa_class_free(class_t *c) 
{
    if (! c) return;
    free(c->class);
    free(c->residue);
    free(c);
}

static result_t* freesasa_result_new(size_t n_atoms)
{
    result_t *r = (result_t*)malloc(sizeof(result_t));

    int nc = freesasa_classify_nclasses();
    int nr = freesasa_classify_nresiduetypes();
    
    r->sasa = (double*)malloc(sizeof(double)*n_atoms);
    r->class = (double*)malloc(sizeof(double)*nc);
    r->residue = (double*)malloc(sizeof(double)*nr);

    for (int i = 0; i < nc; ++i) r->class[i] = 0;
    for (int i = 0; i < nr; ++i) r->residue[i] = 0;
    
    return r;
}

static void freesasa_result_free(result_t *r)
{
    if (! r) return;
    free(r->sasa);
    free(r->class);
    free(r->residue);
    free(r);
}

static void freesasa_get_class_result(freesasa_t *s, freesasa_structure_t *p)
{
    if (s->class) freesasa_class_free(s->class);
    s->class = freesasa_classify_structure(p);
    
    for (size_t i = 0; i < freesasa_structure_n(p); ++i) {
        s->result->class[s->class->class[i]] += s->result->sasa[i];
        s->result->residue[s->class->residue[i]] += s->result->sasa[i];
    }
}

static int freesasa_calc(freesasa_t *s, 
			const freesasa_coord_t *c, const double *r)
{
    struct timeval t1, t2;
    int res = 0;

    gettimeofday(&t1,NULL);
    
    if (s->result) freesasa_result_free(s->result);
    result_t *result;
    s->result = result = freesasa_result_new(s->n_atoms);
    
    switch(s->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        res = freesasa_shrake_rupley(result->sasa, c, r, 
				    s->probe_radius,
				    s->n_sr, s->n_threads);
        break;
    case FREESASA_LEE_RICHARDS:
        res = freesasa_lee_richards(result->sasa, c, r, 
				   s->probe_radius,
				   s->d_lr, s->n_threads);
	break;
    default:
        return freesasa_fail("no SASA algorithm specified.");
    }
    result->total = 0;
    for (int i = 0; i < freesasa_coord_n(c); ++i) {
	result->total += result->sasa[i];
    }

    gettimeofday(&t2,NULL);
    
    s->elapsed_time = (t2.tv_sec-t1.tv_sec); 
    s->elapsed_time += (t2.tv_usec-t1.tv_usec) / 1e6; // s
    
    s->calculated = 1;

    return res;
}

freesasa_t* freesasa_init()
{
    freesasa_t *s = (freesasa_t*) malloc(sizeof(freesasa_t));
    *s = freesasa_def_param;
    return s;
}

void freesasa_copy_param(freesasa_t *target, const freesasa_t *source)
{
    target->alg = source->alg;
    target->n_sr = source->n_sr;
    target->d_lr = source->d_lr;
    target->n_threads = source->n_threads;
    target->probe_radius = source->probe_radius;
    target->calculated = 0;
}

void freesasa_free(freesasa_t *s)
{
    if (! s) return;
    if (s->owns_r) {
        free(s->r);
    }
    if (s->class) freesasa_class_free(s->class);
    if (s->result) freesasa_result_free(s->result);
    if (s->coord) freesasa_coord_free(s->coord);
    if (s->structure) freesasa_structure_free(s->structure);
    free(s);
}

int freesasa_calc_coord(freesasa_t *s, const double *coord, 
                       const double *r, size_t n)
{
    s->calculated = 0;
    s->n_atoms = n;

    // We don't want to store the supplied parameters (to allow const-ness),
    // and want to make sure user doesn't access outdated parameters
    if (s->owns_r && s->r) free(s->r);
    s->owns_r = 0;
    s->r = NULL;
    if (s->coord) free(s->coord);
    s->coord = NULL;

    freesasa_coord_t *c = freesasa_coord_new_linked(coord,n);
    int res = freesasa_calc(s,c,r);
    freesasa_coord_free(c);
    
    return res;
}

int freesasa_calc_pdb(freesasa_t *s, FILE *pdb_file)
{
    s->calculated = 0;
    freesasa_structure_t *p = freesasa_structure_init_from_pdb(pdb_file);
    if (!p) {
        return freesasa_fail("Failure reading PDB-file.");
    }
    s->n_atoms = freesasa_structure_n(p);
    
    if (s->structure) freesasa_structure_free(s->structure);
    s->structure = p;

    //calc OONS radii
    if (s->r) free(s->r);
    s->r = (double*) malloc(sizeof(double)*freesasa_structure_n(p));
    s->owns_r = 1;
    freesasa_structure_r_def(s->r,p);

    int res = freesasa_calc(s,freesasa_structure_xyz(p),s->r);
    if (!res) freesasa_get_class_result(s,p);

    return res;
}
double freesasa_radius(const char* residue_name, const char* atom_name)
{
    return freesasa_classify_radius(residue_name,atom_name);
}
int freesasa_link_coord(freesasa_t *s, const double *coord,
                       double *r, size_t n) 
{
    s->calculated = 0;
    if (s->coord) freesasa_coord_free(s->coord);
    s->coord = freesasa_coord_new_linked(coord,n);

    if (s->r && s->owns_r) free(s->r);
    s->r = r;
    s->owns_r = 0;

    if (s->result) {
	freesasa_result_free(s->result);
	s->result = NULL;
    }

    s->n_atoms = n;

    return FREESASA_SUCCESS;
}

int freesasa_refresh(freesasa_t *s)
{
    s->calculated = 0;
    if (! s->coord ) 
        return freesasa_fail("freesasa_refresh(1) called, but no coordinates available.");
    if (! s->r )
        return freesasa_fail("freesasa_refresh(1) called, but no atomic radii specified.");
    return freesasa_calc(s,s->coord,s->r);
}

int freesasa_set_algorithm(freesasa_t *s, freesasa_algorithm alg)
{
    s->calculated = 0;
    if (alg == FREESASA_SHRAKE_RUPLEY || alg == FREESASA_LEE_RICHARDS) { 
        s->alg = alg;
        return FREESASA_SUCCESS;
    }
    return freesasa_warn("undefined algorithm selected, proceeding with previously selected algorithm (or default).");
}

freesasa_algorithm freesasa_get_algorithm(const freesasa_t *s)
{
    return s->alg;
}

const char* freesasa_algorithm_name(const freesasa_t *s)
{
    assert(s->alg == FREESASA_SHRAKE_RUPLEY || s->alg == FREESASA_LEE_RICHARDS);
    return freesasa_alg_names[s->alg];
}

int freesasa_set_probe_radius(freesasa_t *s,double r) 
{
    s->calculated = 0;
    if (r < 0 || !isfinite(r)) {
        return freesasa_warn("Probe radius r = %f not allowed.",r);
    }
    s->probe_radius = r;
    return FREESASA_SUCCESS;
}

double freesasa_get_probe_radius(const freesasa_t *s)
{
    return s->probe_radius;
}

int freesasa_set_sr_points(freesasa_t *s, int n) {
    s->calculated = 0;
    if (freesasa_srp_n_is_valid(n)) { 
        s->n_sr = n;
        return FREESASA_SUCCESS;
    }
    freesasa_warn("number of test-points must be one of the following:  ");
    if (freesasa_verbosity == FREESASA_V_NORMAL) {
	freesasa_srp_print_n_opt(stderr);
    }
    freesasa_warn("proceeding with default number of test-points (%d)",FREESASA_DEF_SR_N);
    s->n_sr = FREESASA_DEF_SR_N;
    return FREESASA_WARN;
}

int freesasa_get_sr_points(const freesasa_t* s)
{
    if (s->alg == FREESASA_SHRAKE_RUPLEY) return s->n_sr;
    return FREESASA_WARN;
}

int freesasa_set_lr_delta(freesasa_t *s, double d)
{
    s->calculated = 0;
    if (d > 0 && isfinite(d)) {
        s->d_lr = d;
        if (d <= 2) return FREESASA_SUCCESS;
        if (d > 2) {
            freesasa_warn("for regular SASA calculations, " 
                         "slice width should be less than 2 A.");
            freesasa_warn("proceeding with selected value (%f), "
                         "but results might be inaccurate.", d);
            return FREESASA_WARN;
        }
    }
    freesasa_warn("slice width %f invalid.",d);
    freesasa_warn("proceeding with default value (%f).",FREESASA_DEF_LR_D);
    s->d_lr = FREESASA_DEF_LR_D;
    return FREESASA_WARN;
}

double freesasa_get_lr_delta(const freesasa_t *s)
{
    if (s->alg == FREESASA_LEE_RICHARDS) return s->d_lr;
    return -1.0;
}

#if HAVE_LIBPTHREAD
int freesasa_set_nthreads(freesasa_t *s,int n)
{
    if ( n <= 0) {
	s->n_threads = DEF_NTHREADS;
        return freesasa_warn("Number of threads has to be positive. "
			   "Proceeding with default number (%d).", DEF_NTHREADS);
    }
    s->n_threads = n;
    return FREESASA_SUCCESS;
}

int freesasa_get_nthreads(const freesasa_t *s)
{
    return s->n_threads;
}
#endif

size_t freesasa_n_atoms(const freesasa_t *s)
{
    return s->n_atoms;
}

double freesasa_area_total(const freesasa_t *s)
{
    if (! s->calculated) {
        freesasa_fail("SASA calculation has not been performed, "
                     "no total SASA value available.");
        return -1.0;
    }
    return s->result->total;
}
double freesasa_area_class(const freesasa_t* s, freesasa_class c)
{    
    if (! s->calculated) {
        freesasa_fail("SASA calculation has not been performed, "
                     "no total SASA value available.");
	return -1.0;
    }
    assert(c >= FREESASA_POLAR && c <= FREESASA_CLASS_UNKNOWN &&
           "Invalid arguments to freesasa_area_class(2)");
    return s->result->class[c];
}
double freesasa_area_atom(const freesasa_t *s, int i) 
{
    if ( !s->calculated ) {
	freesasa_fail("SASA calculation has not been performed, "
                     "no SASA values available.");
	return -1.0;
    }
    if (i < 0 || i > s->n_atoms-1) {
        freesasa_fail("Atom index %d invalid.",i);
        return -1.0;
    }
    return s->result->sasa[i];
}

const double* freesasa_area_atom_array(const freesasa_t *s)
{
    if ( !s->calculated ) {
        freesasa_fail("SASA calculation has not been performed, "
                     "no SASA values available.");
        return NULL;
    }
    return s->result->sasa;
}
double freesasa_radius_atom(const freesasa_t *s, int i)
{
    if (s->r == NULL) {
        freesasa_fail("No atomic radii have been assigned.");
        return -1.0;
    }
    if (i < 0 || i > s->n_atoms-1) {
        freesasa_fail("Atom index %d invalid.",i);
        return -1.0;
    }
    return s->r[i];
}
const double* freesasa_radius_atom_array(const freesasa_t *s)
{
    if (s->r == NULL) {
        freesasa_fail("No atomic radii have been assigned.");
        return NULL;
    }
    return s->r;
}
void freesasa_set_proteinname(freesasa_t *s,const char *name)
{
    int n;
    if ((n = strlen(name)) > FREESASA_NAME_LIMIT) {
        strcpy(s->proteinname,"...");
        sprintf(s->proteinname+3,"%.*s",FREESASA_NAME_LIMIT-3,
                name+n+3-FREESASA_NAME_LIMIT);
    } else {
        strcpy(s->proteinname,name);
    }
}

const char* freesasa_get_proteinname(const freesasa_t *s) 
{
    return s->proteinname;
}

int freesasa_log(FILE *log, const freesasa_t *s)
{
    if (! s->calculated) {
        const char *msg = "freesasa_log(2) called, but no calculation "
            "has been performed.";
        freesasa_warn(msg);
	if (freesasa_verbosity == FREESASA_V_NORMAL) {
	    fprintf(log,"%s\n",msg);
	}
        return FREESASA_WARN;
    }
    fprintf(log,"name: %s\n",s->proteinname);
    fprintf(log,"algorithm: %s\n",freesasa_alg_names[s->alg]);
    fprintf(log,"probe-radius: %f A\n",s->probe_radius);
#if HAVE_LIBPTHREAD
    fprintf(log,"n_thread: %d\n",s->n_threads);
#endif
    
    switch(s->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        fprintf(log,"n_testpoint: %d\n",s->n_sr);
        break;
    case FREESASA_LEE_RICHARDS:
        fprintf(log,"d_slice: %f Å\n",s->d_lr);
        break;
    default:
        break;
    }
    fprintf(log,"time_elapsed: %f s\n",s->elapsed_time);
    fprintf(log,"n_atoms: %d\n", s->n_atoms);
    return FREESASA_SUCCESS;
}

int freesasa_per_residue(FILE *output, const freesasa_t *s)
{
    if (! output) {
	return freesasa_fail("freesasa_per_residue(2) output file is "
			    "a NULL pointer.");
    }
    if (! s->calculated) {
        return freesasa_fail("freesasa_per_residue(2) called, "
                            "but no cqalculation has been performed.");
    }
    for (int i = 0; i < freesasa_classify_nresiduetypes(); ++i) {
        double sasa = s->result->residue[i];
        if (i < 20 || sasa > 0) {
            fprintf(output,"%s %8.2f\n",freesasa_classify_residue2str(i),sasa);
        }
    }
    return FREESASA_SUCCESS;
}

double freesasa_area_residue(const freesasa_t *s, const char *res_name)
{
    if (! s->calculated) {
        freesasa_warn("freesasa_area_residue(2) called, "
                     "but no calculation has been performed.");
        return -1.0;
    }
    int res = freesasa_classify_residue(res_name);
    return s->result->residue[res];
}

int freesasa_write_pdb(FILE *output, const freesasa_t *s)
{
    if (!s->calculated) {
        return freesasa_fail("Cannot output results before "
                            "calculation has been performed.");
    }
    if (!s->structure) {
        return freesasa_fail("Cannot write B-factors: "
                            "no input structure.");
    }
    if (freesasa_structure_n(s->structure) != s->n_atoms) {
        return freesasa_fail("Cannot write B-factors: "
                            "number of atoms mismatch.");
    }
    if (!s->result->sasa) {
        return freesasa_fail("Cannot write B-factors: "
                            "no SASA values available.");
    }
    return freesasa_structure_write_pdb_bfactors(output,	s->structure,
                                                s->result->sasa);
}
