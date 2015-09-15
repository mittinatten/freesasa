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

// needed for getline()
#define _GNU_SOURCE

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "structure.h"
#include "pdb.h"
#include "classify.h"
#include "coord.h"
#include "util.h"

struct atom {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char descriptor[FREESASA_STRUCTURE_DESCRIPTOR_STRL];
    char line[PDB_LINE_STRL+1];
    char chain_label;
};

struct freesasa_structure {
    struct atom *a;
    freesasa_coord *xyz;
    int number_atoms;
    int number_residues;
    int number_chains;
    int model; // model number
    char *chains; //all chain labels found (as string)
    int *res_first_atom; // first atom of each residue
    char **res_desc;
};

freesasa_structure*
freesasa_structure_new(void)
{
    freesasa_structure *p = malloc(sizeof(freesasa_structure));
    if (p) {
        p->number_atoms = 0;
        p->number_residues = 0;
        p->number_chains = 0;
        p->model = 0;
        p->chains = NULL;
        p->a = NULL;
        p->xyz = NULL;
        p->res_first_atom = NULL;
        p->res_desc = NULL;
        if ((p->chains = malloc(1)) == NULL ||
            (p->a = malloc(sizeof(struct atom))) == NULL ||
            (p->xyz = freesasa_coord_new()) == NULL ||
            (p->res_first_atom = malloc(sizeof(int))) == NULL ||
            (p->res_desc = malloc(sizeof(char*))) == NULL) {
            freesasa_structure_free(p);
            p = NULL;
        } else {
            p->chains[0] = '\0';
        }
    }
    return p;
}

void
freesasa_structure_free(freesasa_structure *p)
{
    if (p == NULL) return;
    if (p->a) free(p->a);
    if (p->xyz) freesasa_coord_free(p->xyz);
    if (p->res_first_atom) free(p->res_first_atom);
    if (p->chains) free(p->chains);
    if (p->res_desc) {
        for (int i = 0; i < p->number_residues; ++i) {
            if (p->res_desc[i]) free(p->res_desc[i]);
        }
        free(p->res_desc);
    }
    free(p);
}

/* returns alt_label if there is any */
char
freesasa_structure_get_pdb_atom(struct atom *a,
                                double *xyz,
                                const char *line)
{
    assert(a); assert(xyz); assert(line);
    a->chain_label = freesasa_pdb_get_chain_label(line);
    freesasa_pdb_get_coord(xyz, line);
    freesasa_pdb_get_atom_name(a->atom_name, line);
    freesasa_pdb_get_res_name(a->res_name, line);
    freesasa_pdb_get_res_number(a->res_number, line);
    return freesasa_pdb_get_alt_coord_label(line);
}

/**
    Handles the reading of pdb-files, returns NULL if problems reading
    or input or malloc failure. Error-messages should explain what
    went wrong.
 */
static freesasa_structure*
from_pdb_impl(FILE *pdb_file,
              struct file_interval it,
              int options)
{
    assert(pdb_file);
    freesasa_structure *p = freesasa_structure_new();
    if (p == NULL) { mem_fail(); return NULL; }
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    char the_alt = ' ';
    fseek(pdb_file,it.begin,SEEK_SET);
    
    while (getline(&line, &len, pdb_file) != -1 && ftell(pdb_file) <= it.end) {
        if (strncmp("ATOM",line,4)==0 ||
            ( (options & FREESASA_INCLUDE_HETATM) &&
              (strncmp("HETATM",line,6) == 0) )
            ) {
            double v[3];
            struct atom a;
            char alt;
            if (freesasa_pdb_ishydrogen(line) &&
                !(options & FREESASA_INCLUDE_HYDROGEN)) continue;
            alt = freesasa_structure_get_pdb_atom(&a,v,line);
            if ((alt != ' ' && the_alt == ' ') || (alt == ' ')) the_alt = alt;
            else if (alt != ' ' && alt != the_alt)              continue;
            if (freesasa_structure_add_atom(p,a.atom_name,a.res_name,a.res_number,
                                            a.chain_label, v[0], v[1], v[2])
                != FREESASA_SUCCESS) { 
                freesasa_structure_free(p);
                p = NULL; 
                break;
            }
            strncpy(p->a[p->number_atoms-1].line,line,PDB_LINE_STRL);
        }

        if (! (options & FREESASA_JOIN_MODELS)) {
            if (strncmp("MODEL",line,5)==0) {
                sscanf(line+10,"%d",&p->model);
            }
            if (strncmp("ENDMDL",line,6)==0) {
                if (p->number_atoms == 0) {
                    freesasa_fail("input had ENDMDL before first ATOM entry.");
                    freesasa_structure_free(p);
                    p = NULL;
                }
                break;
            }
        }
    }
    free(line);
    if (p != NULL && p->number_atoms == 0) {
        freesasa_fail("input had no ATOM entries.");
        freesasa_structure_free(p);
        return NULL;
    }
    return p;
}

static struct file_interval
get_whole_file(FILE* file)
{
    assert(file != NULL);
    struct file_interval interval;
    rewind(file);
    interval.begin = ftell(file);
    fseek(file,0,SEEK_END);
    interval.end = ftell(file);
    rewind(file);
    assert(interval.begin <= interval.end);
    return interval;
}
/**
    Finds the location of all MODEL entries in the file pdb, returns the number
    of models found. *intervals points to a dynamically allocated file intervals
    for each model.

    Returns 0 if no MODEL lines were found (for example an X-ray structure with
    only one model) and sets *intervals to NULL. That means that a return value
    of 0 doesn't have to mean the file is empty.

    Return FREESASA_FAIL if malloc-failure.
*/
static int
get_models(FILE *pdb,
           struct file_interval** intervals)
{
    assert(pdb != NULL);
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    int n = 0, n_end = 0, error = 0;
    long last_pos = ftell(pdb);
    struct file_interval *it = NULL;

    while (getline(&line, &len, pdb) != -1) {
        if (strncmp("MODEL",line,5)==0) {
            ++n;
            it = realloc(it,sizeof(struct file_interval)*n);
            if (!it) return mem_fail();
            it[n-1].begin = last_pos;
        }
        if (strncmp("ENDMDL",line,6)==0) {
            ++n_end;
            if (n != n_end) {
                error = freesasa_fail("%s: Mismatch between MODEL and ENDMDL "
                                      "in input\n",__func__);
                break;
            }
            it[n-1].end = ftell(pdb);
        }
        last_pos = ftell(pdb);
    }
    if (n == 0) { // when there are no models, the whole file is the model
        free(it);
        it = NULL;
    }
    if (error == FREESASA_FAIL) {
        free(it);
        intervals = NULL;
        return FREESASA_FAIL;
    }
    *intervals = it;
    return n;
}

static int
get_chains(FILE *pdb,
           struct file_interval model,
           struct file_interval **intervals,
           int options)
{
    assert(pdb);
    assert(intervals);
    // it is assumed that 'model' is valid for 'pdb'

    int n_chains = 0;
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    struct file_interval *chains = NULL;
    char last_chain = '\0';
    long last_pos = model.begin;
    *intervals = NULL;

    // for each model, find file intervals for each chain, store them
    // in the dynamically growing array chains
    fseek(pdb,model.begin,SEEK_SET);
    while (getline(&line, &len, pdb) != -1 &&
           ftell(pdb) <= model.end ) {
        if (strncmp("ATOM",line,4)==0 ||( (options & FREESASA_INCLUDE_HETATM) &&
                                          (strncmp("HETATM",line,6) == 0) ) ) {
            char chain = freesasa_pdb_get_chain_label(line);
            if (chain != last_chain) {
                if (n_chains > 0) chains[n_chains-1].end = last_pos;
                ++n_chains;
                chains = realloc(chains,sizeof(struct file_interval)*n_chains);
                if (!chains) return mem_fail();
                chains[n_chains-1].begin = last_pos;
                last_chain = chain;
            }
        }
        last_pos = ftell(pdb);
    }
    if (n_chains > 0) {
        chains[n_chains-1].end = last_pos;
        chains[0].begin = model.begin; //preserve model info
        *intervals = chains;
    } else {
        *intervals = NULL;
    }
    return n_chains;
}

freesasa_structure*
freesasa_structure_from_pdb(FILE *pdb_file,
                            int options)
{
    assert(pdb_file);
    return from_pdb_impl(pdb_file, get_whole_file(pdb_file), options);
}

freesasa_structure**
freesasa_structure_array(FILE *pdb,
                         int *n,
                         int options)
{
    assert( (options & FREESASA_SEPARATE_MODELS) ||
            (options & FREESASA_SEPARATE_CHAINS) );
    assert(pdb);
    assert(n);

    struct file_interval *models = NULL;
    struct file_interval whole_file = get_whole_file(pdb);
    int n_models = get_models(pdb,&models), n_chains = 0;
    freesasa_structure **ss = NULL;
    int err = 0;

    if (n_models == FREESASA_FAIL) {
        freesasa_fail("%s: problems reading PDB-file.",__func__);
        return NULL;
    } else if (n_models == 0) {
        models = &whole_file;
        n_models = 1;
    }

    //only keep first model if option not provided
    if (! (options & FREESASA_SEPARATE_MODELS) ) n_models = 1;

    //for each model read chains if requested
    if (options & FREESASA_SEPARATE_CHAINS) {
        for (int i = 0; i < n_models; ++i) {
            struct file_interval* chains = NULL;
            int new_chains = get_chains(pdb,models[i],&chains,options);
            if (new_chains == FREESASA_FAIL) { mem_fail(); ++err; break; }
            if (new_chains == 0) {
                freesasa_warn("%s: warning: No chains found (in model %d).",__func__,i+1);
                continue;
            }
            ss = realloc(ss,sizeof(freesasa_structure*)*(n_chains+new_chains));
            if (!ss) { mem_fail(); return NULL; } // now way of cleaning up, might as well exit
            for (int j = 0; j < new_chains; ++j) {
                freesasa_structure *s = from_pdb_impl(pdb,chains[j],options);
                if (s != NULL) {
                    ss[n_chains+j] = s;
                    ss[n_chains+j]->model = ss[n_chains]->model; // all have the same model number
                } else { ++err; break; }
            }
            n_chains += new_chains;
            free(chains);
        }
        *n = n_chains;
    } else {
        ss = malloc(sizeof(freesasa_structure*)*n_models);
        if (!ss) { mem_fail(); return NULL; }
        for (int i = 0; i < n_models; ++i) {
            freesasa_structure *s = from_pdb_impl(pdb,models[i],options);
            if (s != NULL) ss[i] = s;
            else {
                ++err; 
                break;
            }
        }
        *n = n_models;
    }
    if (err || *n == 0) {
        *n = 0;
        free(ss);
        freesasa_fail("%s: error: Problems reading input.",__func__);
        ss = NULL;
    }
    if (models != &whole_file) free(models);
    return ss;
}

static int
structure_alloc_one(freesasa_structure *p)
{
    assert(p);
    int na = ++p->number_atoms;
    p->a = realloc(p->a,sizeof(struct atom)*na);
    if (p->a != NULL) return FREESASA_SUCCESS;
    return mem_fail();
}

static int
structure_add_chain(freesasa_structure *p, char chain_label)
{
    if (strchr(p->chains,chain_label) == NULL) {
        int n = ++p->number_chains;
        p->chains = realloc(p->chains,n + 1);
        if (p->chains) {
            p->chains[n-1] = chain_label;
            p->chains[n] = '\0';
                assert (strlen(p->chains) == p->number_chains);
        } else {
            return mem_fail();
        }
    }
    return FREESASA_SUCCESS;
}

static inline void
structure_set_atom(struct atom* a, 
                   const char *atom_name,
                   const char *residue_name,
                   const char *residue_number,
                   char chain_label)
{
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    sprintf(a->descriptor,"%c %s %s %s",chain_label,residue_number,residue_name,atom_name);
    a->chain_label = chain_label;
    a->line[0] = '\0';
}

int
freesasa_structure_add_atom(freesasa_structure *p,
                            const char *atom_name,
                            const char *residue_name,
                            const char *residue_number,
                            char chain_label,
                            double x, double y, double z)
{
    assert(p);
    assert(atom_name); assert(residue_name); assert(residue_number);

    int na;

    // check input for consistency
    if (freesasa_classify_validate_atom(residue_name,atom_name) != FREESASA_SUCCESS)
        return freesasa_warn("%s: Skipping atom '%s' in residue '%s'",
                             __func__,atom_name,residue_name);

    // allocate memory, increase number of atoms counter, add chain
    if (structure_alloc_one(p)) return mem_fail();
    if (freesasa_coord_append_xyz(p->xyz, &x, &y, &z, 1)) return mem_fail();
    if (structure_add_chain(p, chain_label)) return mem_fail();

    na = p->number_atoms;
    assert(na > 0);
    structure_set_atom(p->a + na-1,atom_name,residue_name,residue_number,chain_label);

    /* here we assume atoms are ordered sequentially, i.e. residues are
       not mixed in input: if two sequential atoms have different
       residue numbers, a new residue is assumed to begin */
    if (p->number_residues == 0) {
        ++p->number_residues;
        p->res_first_atom[0] = 0;
        if (!(p->res_desc[0] = malloc(FREESASA_STRUCTURE_DESCRIPTOR_STRL)))
            return mem_fail();
        sprintf(p->res_desc[0],"%c %s %s", chain_label,residue_number,residue_name);
    }
    if (na > 1 && strcmp(residue_number,p->a[na-2].res_number)) {
        int naa = p->number_residues+1;
        if (!(p->res_first_atom = realloc(p->res_first_atom, sizeof(int)*naa)))
            return mem_fail();
        p->res_first_atom[naa-1] = na-1;
        if (!(p->res_desc = realloc(p->res_desc, sizeof(char*)*naa)))
            return mem_fail(); 
        if (!(p->res_desc[naa-1] = malloc(FREESASA_STRUCTURE_DESCRIPTOR_STRL)))
            return mem_fail();
        sprintf(p->res_desc[naa-1], "%c %s %s", chain_label, residue_number, residue_name);
        ++p->number_residues;
    }
    return FREESASA_SUCCESS;
}

freesasa_structure*
freesasa_structure_get_chains(const freesasa_structure *p,
                              const char* chains)
{
    assert(p);
    if (strlen(chains) == 0) return NULL;

    freesasa_structure *new_p = freesasa_structure_new();
    if (!new_p) { mem_fail(); return NULL; }
    new_p->model = p->model;

    for (int i = 0; i < p->number_atoms; ++i) {
        struct atom *ai = &(p->a[i]);
        char c = ai->chain_label;
        if (strchr(chains,c) != NULL) {
            const double *v = freesasa_coord_i(p->xyz,i);
            int res = freesasa_structure_add_atom(new_p, ai->atom_name,
                                                  ai->res_name, ai->res_number,
                                                  c, v[0], v[1], v[2]);
            if (res != FREESASA_SUCCESS) { mem_fail(); return NULL; }
        }
    }
    if (new_p->number_atoms == 0) {
        freesasa_structure_free(new_p);
        new_p = NULL;
    }
    return new_p;
}

const char*
freesasa_structure_chain_labels(const freesasa_structure *structure)
{
    assert(structure);
    return structure->chains;
}

const freesasa_coord*
freesasa_structure_xyz(const freesasa_structure *p)
{
    assert(p);
    return p->xyz;
}

int
freesasa_structure_n(const freesasa_structure *p)
{
    assert(p);
    return p->number_atoms;
}

int
freesasa_structure_n_residues(const freesasa_structure *p)
{
    assert(p);
    return p->number_residues;
}

const char*
freesasa_structure_atom_name(const freesasa_structure *p,
                             int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].atom_name;
}

const char*
freesasa_structure_atom_res_name(const freesasa_structure *p,
                                 int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].res_name;
}

const char*
freesasa_structure_atom_res_number(const freesasa_structure *p,
                                   int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].res_number;
}

char
freesasa_structure_atom_chain(const freesasa_structure *p,
                              int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].chain_label;
}

const char*
freesasa_structure_atom_descriptor(const freesasa_structure *p,
                                   int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].descriptor;
}

int
freesasa_structure_residue_atoms(const freesasa_structure *s,
                                 int r_i,
                                 int *first,
                                 int *last)
{
    assert(s); assert(first); assert(last);
    const int naa = s->number_residues;
    assert(r_i >= 0 && r_i < naa);

    *first = s->res_first_atom[r_i];
    if (r_i == naa-1) *last = s->number_atoms-1;
    else *last = s->res_first_atom[r_i+1]-1;
    assert(*last >= *first);
    return FREESASA_SUCCESS;
}

const char*
freesasa_structure_residue_descriptor(const freesasa_structure *s,
                                                  int r_i)
{
    assert(s);
    assert(r_i < s->number_residues);
    return s->res_desc[r_i];
}

int
freesasa_write_pdb(FILE *output,
                       freesasa_result *result,
                       const freesasa_structure *p,
                       const double *radii)
{
    assert(p);
    assert(output);
    assert(radii);
    assert(result);
    assert(result->sasa);

    const double* values = result->sasa;
    char buf[PDB_LINE_STRL+1], buf2[6];
    int n = freesasa_structure_n(p);
    if (p->model > 0) fprintf(output,"MODEL     %4d\n",p->model);
    else fprintf(output,             "MODEL        1\n");
    // Write ATOM entries
    for (int i = 0; i < n; ++i) {
        if (p->a[i].line[0] == '\0') {
            return freesasa_fail("%s: PDB input not valid or not present.",
                                 __func__);
        }
        strncpy(buf,p->a[i].line,PDB_LINE_STRL);
        sprintf(&buf[54],"%6.2f%6.2f",radii[i],values[i]);
        errno = 0;
        if (fprintf(output,"%s\n",buf) < 0)
            return freesasa_fail("%s: %s", __func__, strerror(errno));
    }
    // Write TER  and ENDMDL lines
    errno = 0;
    strncpy(buf2,&buf[6],5);
    buf2[5]='\0';
    if (fprintf(output,"TER   %5d     %4s %c%4s\nENDMDL\n",
                atoi(buf2)+1, p->a[n-1].res_name,
                p->a[n-1].chain_label, p->a[n-1].res_number) < 0) {
        freesasa_fail("%s: %s", __func__, strerror(errno));
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;

}

int
freesasa_structure_model(const freesasa_structure *structure)
{
    return structure->model;
}
