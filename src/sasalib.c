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

#include "sasa.h"
#include "pdbutil.h"
#include "sasalib.h"
#include "srp.h"
#include "oons.h"

#define NBUF 100

struct sasalib_t {
    sasalib_algorithm alg;
    double *sasa;
    double *r;
    protein *p;
    int n_sr;
    double d_lr;
    int n_threads;
    double elapsed_time;
    int owns_protein;
    int customized;
    int inited;
    double total;
    double polar;
    double apolar;
    char proteinname[SASALIB_NAME_LIMIT];
}; 

const sasalib_t sasalib_def_param = {
    .alg = SHRAKE_RUPLEY,
    .sasa = NULL,
    .r = NULL,
    .p = NULL,
    .n_sr = SASALIB_DEF_SR_N,
    .d_lr = SASALIB_DEF_LR_D,
    .n_threads = 1,
    .elapsed_time = 0,
    .owns_protein = 0,
    .customized = 0,
    .inited = 0,
    .total = 0.,
    .polar = 0.,
    .apolar = 0.,
    .proteinname = "undef"
};

const char *sasalib_alg_names[] = {"Lee & Richards", "Shrake & Rupley"};

void sasalib_get_classes(sasalib_t *s)
{
    atomclassifier ac = oons_classes();
    protein *p = s->p;
    int nc = ac.nclasses;
    double sasa_class[nc];
    s->total = 0.;
    for (int i = 0; i < nc; ++i) {
        sasa_class[i] = 0;
    }
    for (size_t i = 0; i < protein_n(p); ++i) {
        int class = ac.classify(protein_atom_res_name(p,i),
                                protein_atom_name(p,i));
        sasa_class[class] += s->sasa[i];
        s->total += s->sasa[i];
    }
    s->polar = sasa_class[polar];
    s->apolar = sasa_class[apolar];
}

int sasalib_calc(sasalib_t *s)
{
    struct timeval t1, t2;
    int res;

    s->r = (double*) malloc(sizeof(double)*protein_n(s->p));
    s->sasa = (double*) malloc(sizeof(double)*protein_n(s->p));
    //calc OONS radii
    protein_r_def(s->r,s->p);

    gettimeofday(&t1,NULL);
    
    switch(s->alg) {
    case SHRAKE_RUPLEY:
        res = sasa_shrake_rupley(s->sasa,protein_xyz(s->p),
                                 s->r,protein_n(s->p),
                                 s->n_sr,s->n_threads);
        break;
    case LEE_RICHARDS:
        res = sasa_lee_richards(s->sasa,protein_xyz(s->p),
                                s->r,protein_n(s->p),
                                s->d_lr,s->n_threads);
        break;
    default:
        fprintf(stderr,"Error: no SASA algorithm specified.\n");
        return 1;
    }
    
    sasalib_get_classes(s);
    
    gettimeofday(&t2,NULL);
    
    s->elapsed_time = (t2.tv_sec-t1.tv_sec); 
    s->elapsed_time += (t2.tv_usec-t1.tv_usec) / 1e6; // s
    
    return res;
}

sasalib_t* sasalib_init()
{
    sasalib_t *s = (sasalib_t*) malloc(sizeof(sasalib_t));
    *s = sasalib_def_param;
    s->inited = 1;
    return s;
}

void sasalib_copy_param(sasalib_t *target, const sasalib_t *source)
{
    target->alg = source->alg;
    target->n_sr = source->n_sr;
    target->d_lr = source->d_lr;
    target->n_threads = source->n_threads;
}

void sasalib_free(sasalib_t *s)
{
    if (s->owns_protein) {
        protein_free(s->p);
    }
    free(s->sasa);
    free(s->r);
    if (s->inited) free(s);
}

int sasalib_calc_pdb(sasalib_t *s, FILE *pdb_file)
{
    s->p = protein_init_from_pdb(pdb_file);
    s->owns_protein = 1;
    return sasalib_calc(s);
}

int sasalib_calc_protein(sasalib_t *s, protein* p)
{
    s->p = p;
    s->owns_protein = 0;
    return sasalib_calc(s);
}

int sasalib_set_algorithm(sasalib_t *s, sasalib_algorithm alg)
{
    if (alg == SHRAKE_RUPLEY || alg == LEE_RICHARDS) { 
        s->alg = alg;
        return 0;
    }
    return 1;
}

sasalib_algorithm sasalib_get_algorithm(const sasalib_t *s)
{
    return s->alg;
}

const char* sasalib_algorithm_name(const sasalib_t *s)
{
    assert(s->alg == SHRAKE_RUPLEY || s->alg == LEE_RICHARDS);
    return sasalib_alg_names[s->alg];
}

int sasalib_set_sr_points(sasalib_t *s, int n) {
    if (srp_n_is_valid(n)) { 
        s->n_sr = n;
        return 0;
    }
    fprintf(stderr,"Error: Number of test-points has to be"
	    " one of the following values:\n  ");
    srp_print_n_opt(stderr);
    fprintf(stderr,"       Proceeding with default value (%d)\n",
            SASALIB_DEF_SR_N);
    s->n_sr = SASALIB_DEF_SR_N;
    return 1;
}

int sasalib_get_sr_points(const sasalib_t* s)
{
    if (s->alg == SHRAKE_RUPLEY) return s->n_sr;
    return -1;
}

int sasalib_set_lr_delta(sasalib_t *s, double d)
{
    if (d > 0 && d < 5.01) {
        s->d_lr = d;
        return 0;
    }
    fprintf(stderr,"Error: slice width has to lie between 0 and 5 Å\n");
    fprintf(stderr,"       Proceeding with default value (%f)\n",
            SASALIB_DEF_LR_D);
    s->d_lr = SASALIB_DEF_LR_D;
    return 1;
}

int sasalib_get_lr_delta(const sasalib_t *s)
{
    if (s->alg == LEE_RICHARDS) return s->d_lr;
    return -1.0;
}

#ifdef PTHREADS
int sasalib_set_nthreads(sasalib_t *s,int n)
{
    if ( n <= 0) {
        fprintf(stderr,"Error: Number of threads has to be positive.\n");
        return 1;
    }
    s->n_threads = n;
    return 0;
}

int sasalib_get_nthreads(const sasalib_t *s)
{
    return s->n_threads;
}
#endif

double sasalib_total(const sasalib_t *s)
{
    return s->total;
}

double sasalib_polar(const sasalib_t *s)
{
    return s->polar;
}

double sasalib_apolar(const sasalib_t *s)
{
    return s->apolar;
}

double sasalib_sasa_atom(const sasalib_t *s, int i)
{
    if (i >= 0 && i < protein_n(s->p)) {
        return s->sasa[i];
    }
    fprintf(stderr,"Error: Atom index %d invalid.\n",i);
    return -1.0;
}

double sasalib_radius_atom(const sasalib_t *s, int i)
{
    if (i >= 0 && i < protein_n(s->p)) {
        return s->r[i];
    }
    fprintf(stderr,"Error: Atom index %d invalid.\n",i);
    return -1.0;
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
    fprintf(log,"# Using van der Waals radii and atom classes defined \n"
            "# by Ooi et al (PNAS 1987, 84:3086-3090) and a probe radius\n"
            "# of %f Å.\n\n", SASA_PROBE_RADIUS);
    fprintf(log,"name: %s\n",s->proteinname);
    fprintf(log,"algorithm: %s\n",sasalib_alg_names[s->alg]);
#ifdef PTHREADS
    fprintf(log,"n_thread: %d\n",s->n_threads);
#endif
    
    switch(s->alg) {
    case SHRAKE_RUPLEY:
        fprintf(log,"n_testpoint: %d\n",s->n_sr);
        break;
    case LEE_RICHARDS:
        fprintf(log,"d_slice: %f Å\n",s->d_lr);
        break;
    default:
        fprintf(log,"Error: no SASA algorithm specified.\n");
        return 1;
    }
    fprintf(log,"time_elapsed: %f s\n",s->elapsed_time);
    fprintf(log,"n_atoms: %d\n", protein_n(s->p));
    return 0;
}

void sasalib_per_residue(FILE *output, const sasalib_t *s)
{
    double sasa = 0;
    char buf[NBUF] = "", prev_buf[NBUF] = "";
    for (int i = 0; i < protein_n(s->p); ++i) {
	sprintf(buf,"%c_%d_%s",protein_atom_chain(s->p,i),
		atoi(protein_atom_res_number(s->p,i)),
		protein_atom_res_name(s->p,i));
	sasa += s->sasa[i];
	if (strcmp(buf,prev_buf)) {
	    fprintf(output,"%s %f\n",prev_buf,sasa);
	    strcpy(prev_buf,buf);
	    sasa = 0;
	}
    }
    fprintf(output,"%s %f\n",buf,sasa);
}
