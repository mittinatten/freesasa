/*
  Copyright Simon Mitternacht 2013-2015.

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
#include <stdio.h>
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
static freesasa_verbosity verbosity;

typedef struct {
    int *class;
    int *residue_type;
} freesasa_class_data;

typedef struct {
    double *sasa;
    double total;
    double *class;
    double *residue_type;
} freesasa_result;

struct freesasa {
    double *r;
    freesasa_coord *coord;
    freesasa_class_data *class;
    int n_atoms;
    freesasa_structure *structure;

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
    freesasa_result *result;
};

enum error_type {NO_RESULT};

const freesasa freesasa_def_param = {
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
    fflush(stderr);
}

int freesasa_fail(const char *format,...)
{
    if (verbosity == FREESASA_V_SILENT) return FREESASA_FAIL;
    va_list arg;
    va_start(arg, format);
    freesasa_err_impl(FREESASA_FAIL,format,arg);
    va_end(arg);
    return FREESASA_FAIL;
}

int freesasa_warn(const char *format,...)
{
    if (verbosity == FREESASA_V_SILENT) return FREESASA_WARN;
    va_list arg;
    va_start(arg, format);
    freesasa_err_impl(FREESASA_WARN,format,arg);
    va_end(arg);
    return FREESASA_WARN;
}

static freesasa_class_data* freesasa_classify_structure(const freesasa_structure *p)
{
    assert(p);
    assert(freesasa_structure_n(p) >= 0 &&
           "Trying to classify atoms of illegal structure.");

    size_t n = freesasa_structure_n(p);

    freesasa_class_data* c = (freesasa_class_data*)malloc(sizeof(freesasa_class_data));
    assert(c);
    c->class = (int*)malloc(sizeof(int)*n);
    c->residue_type = (int*)malloc(sizeof(int)*n);
    assert(c->class && c->residue_type);
    
    for (int i = 0; i < n; ++i) {
        const char *res_name = freesasa_structure_atom_res_name(p,i);
        const char *atom_name = freesasa_structure_atom_name(p,i);
        c->class[i] = freesasa_classify_class(res_name,atom_name);
        c->residue_type[i] = freesasa_classify_residue(res_name);
    }

    return c;
}

static void freesasa_class_free(freesasa_class_data *c)
{
    if (! c) return;
    free(c->class);
    free(c->residue_type);
    free(c);
}

static freesasa_result* freesasa_result_new(size_t n_atoms)
{
    freesasa_result *r = (freesasa_result*)malloc(sizeof(freesasa_result));
    assert(r);
    int nc = freesasa_classify_nclasses();
    int nr = freesasa_classify_nresiduetypes();

    r->sasa = (double*)malloc(sizeof(double)*n_atoms);
    r->class = (double*)malloc(sizeof(double)*nc);
    r->residue_type = (double*)malloc(sizeof(double)*nr);
    assert(r->sasa && r->class && r->residue_type);
    
    for (int i = 0; i < nc; ++i) r->class[i] = 0;
    for (int i = 0; i < nr; ++i) r->residue_type[i] = 0;

    return r;
}

static void freesasa_result_free(freesasa_result *r)
{
    if (! r) return;
    free(r->sasa);
    free(r->class);
    free(r->residue_type);
    free(r);
}

static void freesasa_get_class_result(freesasa *s, freesasa_structure *p)
{
    assert(s);
    assert(p);
    
    if (s->class) freesasa_class_free(s->class);
    s->class = freesasa_classify_structure(p);

    for (size_t i = 0; i < freesasa_structure_n(p); ++i) {
        s->result->class[s->class->class[i]] += s->result->sasa[i];
        s->result->residue_type[s->class->residue_type[i]] += s->result->sasa[i];
    }
}

static int freesasa_calc(freesasa *s,
                         const freesasa_coord *c, const double *r)
{
    assert(s);
    assert(c);
    assert(r);
    
    struct timeval t1, t2;
    int res = 0;

    gettimeofday(&t1,NULL);

    if (s->result) freesasa_result_free(s->result);
    freesasa_result *result;
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
        return freesasa_fail("%s: no SASA algorithm specified.", __func__);
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

freesasa* freesasa_new()
{
    freesasa *s = (freesasa*) malloc(sizeof(freesasa));
    assert(s);
    *s = freesasa_def_param;
    return s;
}

void freesasa_copy_param(freesasa *target, const freesasa *source)
{
    assert(target);
    assert(source);
    target->alg = source->alg;
    target->n_sr = source->n_sr;
    target->d_lr = source->d_lr;
    target->n_threads = source->n_threads;
    target->probe_radius = source->probe_radius;
    target->calculated = 0;
}

void freesasa_free(freesasa *s)
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
int freesasa_calc_coord(freesasa *s, const double *coord,
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

    freesasa_coord *c = freesasa_coord_new_linked(coord,n);
    int res = freesasa_calc(s,c,r);
    freesasa_coord_free(c);

    return res;
}
// assumes s has a structure
static int freesasa_calc_structure(freesasa *s)
{
    assert(s);
    assert(s->structure);

    //calc OONS radii
    if (s->owns_r && s->r) free(s->r);
    s->n_atoms = freesasa_structure_n(s->structure);
    s->r = (double*) malloc(sizeof(double)*s->n_atoms);
    assert(s->r);
    s->owns_r = 1;
    freesasa_structure_r_def(s->r,s->structure);
    
    int res = freesasa_calc(s,freesasa_structure_xyz(s->structure),s->r);
    if (res == FREESASA_SUCCESS) freesasa_get_class_result(s,s->structure);
    return res;
}
int freesasa_calc_pdb(freesasa *s, FILE *pdb_file)
{
    assert(s);
    assert(pdb_file);
    s->calculated = 0;
    freesasa_structure *p = freesasa_structure_from_pdb(pdb_file);
    if (!p) {
        return freesasa_fail("%s: Failure reading PDB-file.", __func__);
    }
    
    if (s->structure) freesasa_structure_free(s->structure);
    s->structure = p;

    return freesasa_calc_structure(s);
}

int freesasa_calc_atoms(freesasa *s, const double *coord, 
                         const char **resnames, 
                         const char **atomnames, size_t n)
{
    assert(s);
    assert(coord);
    assert(resnames);
    assert(atomnames);
    assert(n > 0);
    int status;
    int warn = 0, fail = 0;

    if (s->structure) freesasa_structure_free(s->structure);
    
    freesasa_structure *p = freesasa_structure_new();
    for (size_t i = 0; i < n; ++i) {
        status = freesasa_structure_add_atom(p,atomnames[i],resnames[i],"   1",'A',
                                              coord[i*3],coord[i*3+1],coord[i*3+2]);
        if (status == FREESASA_WARN) ++warn;
        if (status == FREESASA_FAIL) ++fail;
    } 
    s->structure = p;

    status = freesasa_calc_structure(s);

    if (status == FREESASA_FAIL) ++fail;
    if (fail == 0 && warn == 0) return FREESASA_SUCCESS;
    if (fail == 0) return FREESASA_WARN;

    return FREESASA_FAIL;
}

double freesasa_radius(const char* residue_name, const char* atom_name)
{
    return freesasa_classify_radius(residue_name,atom_name);
}
int freesasa_link_coord(freesasa *s, const double *coord,
                        double *r, size_t n)
{
    assert(s);
    assert(coord);
    assert(n > 0);
    
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

int freesasa_refresh(freesasa *s)
{
    assert(s);
    assert(s->coord);
    assert(s->r);
    s->calculated = 0;
    return freesasa_calc(s,s->coord,s->r);
}

int freesasa_set_algorithm(freesasa *s, freesasa_algorithm alg)
{
    assert(s);
    assert(alg == FREESASA_SHRAKE_RUPLEY || alg == FREESASA_LEE_RICHARDS);
    s->calculated = 0;
    s->alg = alg;
    return FREESASA_SUCCESS;
}

freesasa_algorithm freesasa_get_algorithm(const freesasa *s)
{
    assert(s);
    return s->alg;
}

const char* freesasa_algorithm_name(const freesasa *s)
{
    assert(s);
    assert(s->alg == FREESASA_SHRAKE_RUPLEY || s->alg == FREESASA_LEE_RICHARDS);
    return freesasa_alg_names[s->alg];
}

int freesasa_set_probe_radius(freesasa *s,double r)
{
    assert(s);
    s->calculated = 0;
    if (r < 0 || !isfinite(r)) {
        return freesasa_warn("%s: Probe radius r = %f not allowed.",r,__func__);
    }
    s->probe_radius = r;
    return FREESASA_SUCCESS;
}

double freesasa_get_probe_radius(const freesasa *s)
{
    assert(s);
    return s->probe_radius;
}

int freesasa_set_sr_points(freesasa *s, int n) {
    assert(s);
    s->calculated = 0;
    if (freesasa_srp_n_is_valid(n)) {
        s->n_sr = n;
        return FREESASA_SUCCESS;
    }
    freesasa_warn("%s: number of test-points must be one of the following:  ",
                  __func__);
    if (verbosity == FREESASA_V_NORMAL) {
        freesasa_srp_print_n_opt(stderr);
    }
    freesasa_warn("%s: proceeding with default number of test-points (%d)",
                  __func__, FREESASA_DEF_SR_N);
    s->n_sr = FREESASA_DEF_SR_N;
    return FREESASA_WARN;
}

int freesasa_get_sr_points(const freesasa* s)
{
    assert(s);
    if (s->alg == FREESASA_SHRAKE_RUPLEY) return s->n_sr;
    return FREESASA_WARN;
}

int freesasa_set_lr_delta(freesasa *s, double d)
{
    assert(s);
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

double freesasa_get_lr_delta(const freesasa *s)
{
    if (s->alg == FREESASA_LEE_RICHARDS) return s->d_lr;
    return -1.0;
}

#if HAVE_LIBPTHREAD
int freesasa_set_nthreads(freesasa *s,int n)
{
    if ( n <= 0) {
        s->n_threads = DEF_NTHREADS;
        return freesasa_warn("Number of threads has to be positive. "
                             "Proceeding with default number (%d).", DEF_NTHREADS);
    }
    s->n_threads = n;
    return FREESASA_SUCCESS;
}

int freesasa_get_nthreads(const freesasa *s)
{
    return s->n_threads;
}
#endif

size_t freesasa_n_atoms(const freesasa *s)
{
    return s->n_atoms;
}

double freesasa_area_total(const freesasa *s)
{
    assert(s->calculated);
    assert(s->result);
    return s->result->total;
}
double freesasa_area_class(const freesasa* s, freesasa_class c)
{
    assert(s->calculated);
    assert(c >= FREESASA_POLAR && c <= FREESASA_CLASS_UNKNOWN &&
           "Invalid arguments to freesasa_area_class(2)");
    return s->result->class[c];
}
double freesasa_area_atom(const freesasa *s, int i)
{
    assert(s->calculated);
    assert(i >= 0 || i < s->n_atoms);
    assert(s->result);
    return s->result->sasa[i];
}

const double* freesasa_area_atom_array(const freesasa *s)
{
    assert(s->calculated);
    return s->result->sasa;
}
double freesasa_radius_atom(const freesasa *s, int i)
{
    assert(s->r && "no radii available");
    assert(i >= 0 || i < s->n_atoms);
    return s->r[i];
}
const double* freesasa_radius_atom_array(const freesasa *s)
{
    assert(s->r && "no radii available");
    return s->r;
}
void freesasa_set_proteinname(freesasa *s,const char *name)
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

const char* freesasa_get_proteinname(const freesasa *s)
{
    return s->proteinname;
}

int freesasa_log(const freesasa *s, FILE *log)
{
    assert(s->calculated);
    assert(log);
    fprintf(log,"name: %s\n",s->proteinname);
    fprintf(log,"algorithm: %s\n",freesasa_alg_names[s->alg]);
    fprintf(log,"probe-radius: %f A\n",s->probe_radius);
    if(HAVE_LIBPTHREAD)
        fprintf(log,"n_thread: %d\n",s->n_threads);
    
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

int freesasa_per_residue_type(const freesasa *s, FILE *output)
{
    assert(s->calculated);
    for (int i = 0; i < freesasa_classify_nresiduetypes(); ++i) {
        double sasa = s->result->residue_type[i];
        if (i < 20 || sasa > 0) {
            fprintf(output,"RES: %s %9.2f\n",
                    freesasa_classify_residue2str(i),sasa);
        }
    }
    return FREESASA_SUCCESS;
}

static double freesasa_single_residue_sasa(const freesasa *s, int r_i)
{
    assert(s);
    assert(s->calculated);
    assert(s->structure);
    int first, last;
    const freesasa_structure *p = s->structure;
    const double *sasa = s->result->sasa;            
    freesasa_structure_residue_atoms(p,r_i,&first,&last);
    double a = 0;
    for (int j = first; j <= last; ++j) {
        a += sasa[r_i];
    }        
    return a;
}

int freesasa_per_residue(const freesasa *s, FILE *output)
{
    assert(s);
    assert(s->structure);
    assert(s->calculated);
    
    const freesasa_structure *p = s->structure;
    const int naa = freesasa_structure_n_residues(p);
    for (int i = 0; i < naa; ++i) {
        fprintf(output,"SEQ: %s %7.2f\n",
                freesasa_structure_residue_descriptor(p,i),
                freesasa_single_residue_sasa(s,i));
    }
    return FREESASA_SUCCESS;
}

double freesasa_area_residue(const freesasa *s, const char *res_name)
{
    assert(s->calculated);
    int res = freesasa_classify_residue(res_name);
    return s->result->residue_type[res];
}

int freesasa_write_pdb(const freesasa *s, FILE *output)
{
    assert(freesasa_structure_n(s->structure) == s->n_atoms);
    assert(s->calculated);
    assert(s->structure);
    assert(s->result->sasa);
    return freesasa_structure_write_pdb_bfactors(s->structure,output,
                                                 s->result->sasa);
}


static freesasa_strvp* freesasa_alloc_strvp(size_t n)
{
    const int STRL=FREESASA_STRUCTURE_DESCRIPTOR_STRL;
    freesasa_strvp* svp = (freesasa_strvp*) malloc(sizeof(freesasa_strvp));
    assert(svp);
    svp->value = (double*) malloc(sizeof(double)*n);
    svp->string = (char**) malloc(sizeof(char*)*n);
    assert(svp->value && svp->string);

    svp->n = n;
    for (size_t i = 0; i < n; ++i) {
        svp->string[i] = (char*) malloc(sizeof(char)*STRL);
        assert(svp->string[i]);
    }
    return svp;
}

void freesasa_strvp_free(freesasa_strvp *svp)
{
    if (svp->value) free(svp->value);
    if (svp->string) {
        for (size_t i = 0; i < svp->n; ++i) {
            if (svp->string[i]) free(svp->string[i]);
        }
        free(svp->string);
    }
    free(svp);
}
    
freesasa_strvp* freesasa_string_value_pairs(const freesasa *s,freesasa_result_type type)
{
    assert(s);
    assert(s->calculated);

    const freesasa_structure *p = s->structure;
    freesasa_strvp *svp = NULL;
    size_t n;
    switch (type) {
    case FREESASA_ATOMS:
        n = freesasa_n_atoms(s);
        svp = freesasa_alloc_strvp(n);
        for (size_t i = 0; i < n; ++i) {
            svp->value[i] = freesasa_area_atom(s,i);
            strcpy(svp->string[i],freesasa_structure_atom_descriptor(p,i));
        }
        break;
    case FREESASA_RESIDUES:
        n = freesasa_structure_n_residues(p);
        svp = freesasa_alloc_strvp(n);
        for (size_t i = 0; i < n; ++i) {
            svp->value[i] = freesasa_single_residue_sasa(s,i);
            strcpy(svp->string[i],freesasa_structure_residue_descriptor(p,i));
        }
        break;
    default:
        freesasa_fail("Error: Illegal result type in freesasa_string_value_pairs().\n");
    }
    return svp;
}

int freesasa_set_verbosity(freesasa_verbosity s) {
    if (s == FREESASA_V_NORMAL || s == FREESASA_V_SILENT) {
        verbosity = s;
        return FREESASA_SUCCESS;
    }
    return FREESASA_WARN;
}
freesasa_verbosity freesasa_get_verbosity(void) {
    return verbosity;
}
