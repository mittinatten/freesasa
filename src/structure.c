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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <errno.h>
#include "structure.h"
#include "pdb.h"
#include "classifier.h"
#include "coord.h"
#include "util.h"

struct atom {
    char *res_name;
    char *res_number;
    char *atom_name;
    char *symbol;
    char *descriptor;
    char *line;
    char chain_label;
};

static const struct atom empty_atom =
     {NULL, NULL, NULL, NULL, NULL, NULL, '\0'};

struct freesasa_structure {
    struct atom **a;
    coord_t *xyz;
    double *radius;
    int number_atoms;
    int number_residues;
    int number_chains;
    int model; // model number
    char *chains; //all chain labels found (as string)
    int *res_first_atom; // first atom of each residue
    char **res_desc;
};

static const struct freesasa_structure empty_structure = 
    {NULL,NULL,NULL,0,0,0,0,NULL,NULL,NULL};

static int
guess_symbol(char *symbol,
             const char *name);
static void
atom_free(struct atom *a)
{
    if (a == NULL) return;
    free(a->res_name);
    free(a->res_number);
    free(a->atom_name);
    free(a->symbol);
    free(a->descriptor);
    free(a->line);
    free(a);
}

static struct atom *
atom_new(const char *residue_name,
         const char *residue_number,
         const char *atom_name,
         const char *symbol,
         char chain_label)
{
    struct atom *a = malloc(sizeof(struct atom));
    int dlen = 0;
    if (a == NULL) { 
        mem_fail(); 
        return NULL; 
    }

    *a = empty_atom;

    a->line = NULL;
    a->chain_label = chain_label;

    a->res_name = strdup(residue_name);
    a->res_number = strdup(residue_number);
    a->atom_name = strdup(atom_name);
    a->symbol = strdup(symbol);

    dlen = strlen(residue_number) + strlen(residue_name)
        + strlen(atom_name) + 4;
    a->descriptor = malloc(dlen+1);

    if (!a->res_name || !a->res_number || !a->atom_name ||
        !a->symbol || !a->descriptor) {
        mem_fail();
        atom_free(a);
        return NULL;
    }

    sprintf(a->descriptor,"%c %s %s %s",
            chain_label,residue_number,residue_name,atom_name);

    return a;
}

static struct atom *
atom_new_from_line(const char *line,
                   char *alt_label) 
{
    assert(line);
    const int buflen = strlen(line);
    int flag;
    struct atom *a;
    char aname[buflen], rname[buflen], rnumber[buflen], symbol[buflen];

    if (alt_label) *alt_label = freesasa_pdb_get_alt_coord_label(line);

    freesasa_pdb_get_atom_name(aname, line);
    freesasa_pdb_get_res_name(rname, line);
    freesasa_pdb_get_res_number(rnumber, line);

    flag = freesasa_pdb_get_symbol(symbol, line);
    if (flag == FREESASA_FAIL) guess_symbol(symbol,aname);

    a = atom_new(rname, rnumber, aname, symbol, freesasa_pdb_get_chain_label(line));
    
    if (a == NULL) return NULL;

    a->line = strdup(line);
    
    if (a->line != NULL) {
        return a;
    } else {
        mem_fail();
        atom_free(a);
        return NULL;
    }
}

freesasa_structure*
freesasa_structure_new(void)
{
    freesasa_structure *p = malloc(sizeof(freesasa_structure));

    if (p == NULL) {
        mem_fail();
        return NULL;
    }

    *p = empty_structure;

    if ((p->chains = malloc(1)) == NULL ||
        (p->a = malloc(sizeof(struct atom*))) == NULL ||
        (p->xyz = freesasa_coord_new()) == NULL ||
        (p->res_first_atom = malloc(sizeof(int))) == NULL ||
        (p->res_desc = malloc(sizeof(char*))) == NULL) {
        freesasa_structure_free(p);
        mem_fail();
        return NULL;
    }
    p->chains[0] = '\0';

    return p;
}

void
freesasa_structure_free(freesasa_structure *p)
{
    if (p == NULL) return;
    if (p->a) {
        for (int i = 0; i < p->number_atoms; ++i) 
            if (p->a[i]) atom_free(p->a[i]);
        free(p->a);
    }
    if (p->xyz) freesasa_coord_free(p->xyz);
    if (p->res_desc) {
        for (int i = 0; i < p->number_residues; ++i)
            if (p->res_desc[i]) free(p->res_desc[i]);
       free(p->res_desc);
    }
    free(p->radius);
    free(p->res_first_atom);
    free(p->chains);
    free(p);
}

/**
    This function is called when either the symbol field is missing
    from an ATOM record, or when an atom is added directly using
    freesasa_structure_add_atom() or
    freesasa_structure_add_atom_wopt(). The symbol is in turn only
    used if the atom cannot be recognized by the classifier.
*/
static int
guess_symbol(char *symbol,
             const char *name) 
{
    // if the first position is empty, assume that it is a one letter element
    // e.g. " C  "
    if (name[0] == ' ') { 
        symbol[0] = ' ';
        symbol[1] = name[1];
        symbol[2] = '\0';
    } else { 
        // if the string has padding to the right, it's a
        // two-letter element, e.g. "FE  "
        if (name[3] == ' ') {
            strncpy(symbol,name,2);
            symbol[2] = '\0';
        } else { 
            // If it's a four-letter string, it's hard to say,
            // assume only the first letter signifies the element
            symbol[0] = ' ';
            symbol[1] = name[0];
            symbol[2] = '\0';
            return freesasa_warn("Guessing that atom '%s' is symbol '%s'",
                                 name,symbol);
        }
    }
    return FREESASA_SUCCESS;
}
static int
structure_add_chain(freesasa_structure *p,
                    char chain_label)
{
    if (strchr(p->chains,chain_label) == NULL) {
        int n = ++p->number_chains;
        char *pc = p->chains;
        p->chains = realloc(p->chains,n + 1);
        if (p->chains) {
            p->chains[n-1] = chain_label;
            p->chains[n] = '\0';
            assert (strlen(p->chains) == p->number_chains);
        } else {
            p->chains = pc;
            return mem_fail();
        }
    }
    return FREESASA_SUCCESS;
}

/**
    Get the radius of an atom, and fail, warn and/or guess depending
    on the options.
 */
static int
structure_check_atom_radius(double *radius,
                            struct atom *a,
                            const freesasa_classifier* classifier,
                            int options)
{
    *radius = classifier->radius(a->res_name, a->atom_name, classifier);
    if (*radius < 0) {
        if (options & FREESASA_HALT_AT_UNKNOWN) {
            return freesasa_fail("in %s(): atom '%s %s' unknown.",
                                 __func__, a->res_name, a->atom_name);
        } else if (options & FREESASA_SKIP_UNKNOWN) {
            return freesasa_warn("Skipping unknown atom '%s %s'.",
                                 a->res_name, a->atom_name, a->symbol, *radius);
        } else {
            *radius = freesasa_guess_radius(a->symbol);
            if (*radius < 0) {
                *radius = +0.;
                freesasa_warn("Atom '%s %s' unknown and ",
                              "can't guess radius of symbol '%s'. "
                              "Assigning radius 0 A.",
                              a->res_name, a->atom_name, a->symbol);
            } else {
                freesasa_warn("Atom '%s %s' unknown, guessing element is '%s', "
                              "and radius %.3f A.",
                              a->res_name, a->atom_name, a->symbol, *radius);
            }
            // do not return FREESASA_WARN here, because we will keep the atom
        }
    }
    return FREESASA_SUCCESS;
}

/** adds an atom to the structure */
static int
structure_add_atom(freesasa_structure *p,
                   struct atom *a,
                   double *xyz,
                   const freesasa_classifier* classifier,
                   int options)
{
    assert(p); assert(a); assert(xyz);
    int na, ret;
    double r, *pr = p->radius;
    struct atom **pa = p->a;
    
    if (classifier == NULL) {
        classifier = &freesasa_default_classifier;
    }

    // calculate radius and check if we should keep the atom (based on options)
    ret = structure_check_atom_radius(&r, a, classifier, options);
    if (ret == FREESASA_FAIL) return fail_msg("Halting at unknown atom.");
    if (ret == FREESASA_WARN) return FREESASA_WARN;
    assert(r >= 0);

    // if it's a keeper store the radius
    na = p->number_atoms+1;
    if ((p->radius = realloc(p->radius,sizeof(double)*na)) == NULL) {
        p->radius = pr;
        return mem_fail();
    }
    p->radius[na-1] = r;

    // allocate memory for atom, add chain
    if ((p->a = realloc(p->a,sizeof(struct atom*)*na)) == NULL) {
        p->a = pa;
        return mem_fail();
    }

    if (freesasa_coord_append(p->xyz, xyz, 1)) return mem_fail();
    if (structure_add_chain(p, a->chain_label)) return mem_fail();

    /* here we assume atoms are ordered sequentially, i.e. residues are
       not mixed in input: if two sequential atoms have different
       residue numbers, a new residue is assumed to begin */
    if (p->number_residues == 0) {
        ++p->number_residues;
        p->res_first_atom[0] = 0;
        if (!(p->res_desc[0] = malloc(strlen(a->res_number)+
                                      strlen(a->res_name)+4)))
            return mem_fail();
        sprintf(p->res_desc[0], "%c %s %s",
                a->chain_label, a->res_number, a->res_name);
    }

    if (na > 1 && strcmp(a->res_number,p->a[na-2]->res_number)) {
        int naa = p->number_residues+1, *prfa = p->res_first_atom;
        char **prd = p->res_desc;
        
        if (!(p->res_first_atom = realloc(p->res_first_atom, sizeof(int)*naa))) {
            p->res_first_atom = prfa;
            return mem_fail();
        }
        p->res_first_atom[naa-1] = na-1;
        if (!(p->res_desc = realloc(p->res_desc, sizeof(char*)*naa))) {
            p->res_desc = prd;
            return mem_fail();
        }
        if (!(p->res_desc[naa-1] = malloc(strlen(a->res_number)+
                                          strlen(a->res_name)+4)))
            return mem_fail();
        sprintf(p->res_desc[naa-1], "%c %s %s",
                a->chain_label, a->res_number, a->res_name);
        ++p->number_residues;
    }

    // by doing this last, we can free as much memory as possible if anything fails
    p->a[na-1] = a;
    ++p->number_atoms;

    return FREESASA_SUCCESS;
}

/**
    Handles the reading of PDB-files, returns NULL if problems reading
    or input or malloc failure. Error-messages should explain what
    went wrong.
 */
static freesasa_structure*
from_pdb_impl(FILE *pdb_file,
              struct file_range it,
              const freesasa_classifier *classifier,
              int options)
{
    assert(pdb_file);
    int err = 0;
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    char alt, the_alt = ' ';
    double v[3];
    struct atom *a = NULL;
    freesasa_structure *p = freesasa_structure_new();
 
    if (p == NULL) { 
        mem_fail();
        return NULL;
    }
    
    fseek(pdb_file,it.begin,SEEK_SET);
    
    while (getline(&line, &len, pdb_file) != -1 && ftell(pdb_file) <= it.end) {
        
        if (strncmp("ATOM",line,4)==0 || ( (options & FREESASA_INCLUDE_HETATM) &&
                                           (strncmp("HETATM",line,6) == 0) )) {
            
            if (freesasa_pdb_ishydrogen(line) &&
                !(options & FREESASA_INCLUDE_HYDROGEN))
                continue;
 
            a = atom_new_from_line(line,&alt);
            
            if (a == NULL) {
                fail_msg("");
                ++err;
                break;
            }
            
            if ((alt != ' ' && the_alt == ' ') || (alt == ' '))
                the_alt = alt;
            else if (alt != ' ' && alt != the_alt)
                continue;
                        
            if (freesasa_pdb_get_coord(v,line) == FREESASA_FAIL ||
                structure_add_atom(p,a,v,classifier,options) == FREESASA_FAIL) { 
                fail_msg("");
                atom_free(a);
                ++err;
                break;
            }
        }

        if (! (options & FREESASA_JOIN_MODELS)) {
            if (strncmp("MODEL",line,5)==0)  sscanf(line+10,"%d",&p->model);
            if (strncmp("ENDMDL",line,6)==0) break;
        }
    }
    
    free(line);
    if (err == 0 && p->number_atoms == 0) {
        freesasa_fail("Input had no valid ATOM or HETATM lines.");
        ++err;
    }
    if (err) {
        freesasa_structure_free(p);
        p = NULL;
    }
    return p;
}


int
freesasa_structure_add_atom_wopt(freesasa_structure *p,
                                 const char *atom_name,
                                 const char *residue_name,
                                 const char *residue_number,
                                 char chain_label,
                                 double x, double y, double z,
                                 const freesasa_classifier *classifier,
                                 int options)
{
    assert(p);
    assert(atom_name); assert(residue_name); assert(residue_number);

    struct atom *a;
    char symbol[PDB_ATOM_SYMBOL_STRL+1];
    double v[3] = {x,y,z};
    int ret, warn = 0;

    if (guess_symbol(symbol,atom_name) == FREESASA_WARN &&
        options & FREESASA_SKIP_UNKNOWN) ++warn;

    a = atom_new(residue_name,residue_number,atom_name,symbol,chain_label);
    if (a == NULL) return mem_fail();

    ret = structure_add_atom(p,a,v,classifier,options);

    if (ret == FREESASA_FAIL) {
        atom_free(a);
        return ret;
    }

    if (warn) return FREESASA_WARN;

    return ret;
}

int 
freesasa_structure_add_atom(freesasa_structure *p,
                            const char *atom_name,
                            const char *residue_name,
                            const char *residue_number,
                            char chain_label,
                            double x, double y, double z)
{
    return freesasa_structure_add_atom_wopt(p, atom_name, residue_name, residue_number,
                                            chain_label, x, y, z, NULL, 0);
}

freesasa_structure *
freesasa_structure_from_pdb(FILE *pdb_file,
                            const freesasa_classifier* classifier,
                            int options)
{
    assert(pdb_file);
    return from_pdb_impl(pdb_file, freesasa_whole_file(pdb_file),
                         classifier, options);
}

freesasa_structure **
freesasa_structure_array(FILE *pdb,
                         int *n,
                         const freesasa_classifier *classifier,
                         int options)
{
    assert(pdb);
    assert(n);

    if( ! (options & FREESASA_SEPARATE_MODELS ||
           options & FREESASA_SEPARATE_CHAINS) ) {
        fail_msg("Options need to specify at least one of FREESASA_SEPARATE_CHAINS "
                 "and FREESASA_SEPARATE_MODELS.");
        return NULL;
    }

    struct file_range *models = NULL;
    struct file_range whole_file = freesasa_whole_file(pdb);
    int n_models = freesasa_pdb_get_models(pdb,&models), n_chains = 0;
    freesasa_structure **ss = NULL, **ssb;
    int err = 0;

    if (n_models == FREESASA_FAIL) {
        fail_msg("Problems reading PDB-file.");
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
            struct file_range* chains = NULL;
            int new_chains = freesasa_pdb_get_chains(pdb,models[i],&chains,options);
            int j0 = n_chains;
            if (new_chains == FREESASA_FAIL) { mem_fail(); ++err; break; }
            if (new_chains == 0) {
                freesasa_warn("in %s(): No chains found (in model %d).",__func__,i+1);
                continue;
            }
            ssb = ss;
            ss = realloc(ss,sizeof(freesasa_structure*)*(n_chains + new_chains));

            if (!ss) {
                free(chains);
                ss = ssb; // we'll free it below
                mem_fail();
                ++err;
                break;
            } 

            n_chains += new_chains;

            // to facilitate cleanup if reading fails below
            for (int j = 0; j < new_chains; ++j) ss[j0+j] = NULL;
            
            for (int j = 0; j < new_chains; ++j) {
                freesasa_structure *s = from_pdb_impl(pdb,chains[j],classifier,options);
                ss[j0+j] = s;
                if (s == NULL) {
                    ++err;
                    break;
                }
                // all have the same model number
                ss[j0+j]->model = ss[n_chains-new_chains]->model; 
            }

            free(chains);
        }
        *n = n_chains;
    } else {
        ss = malloc(sizeof(freesasa_structure*)*n_models);
        if (!ss) {
            mem_fail();
            ++err;
        } else {
            // to facilitate cleanup if reading fails below
            for (int i = 0; i < n_models; ++i) ss[i] = NULL; 
            *n = n_models;
            
            for (int i = 0; i < n_models; ++i) {
                freesasa_structure *s = from_pdb_impl(pdb,models[i],classifier,options);
                ss[i] = s;
                if (s == NULL) {
                    ++err; 
                    break;
                }
            }
        }
    }

    // clean up
    if (err || *n == 0) {
        if (ss) for (int i = 0; i < *n; ++i) freesasa_structure_free(ss[i]);
        *n = 0;
        free(ss);
        ss = NULL;
    }
    if (models != &whole_file) free(models);

    return ss;
}

freesasa_structure*
freesasa_structure_get_chains(const freesasa_structure *p,
                              const char* chains)
{
    assert(p);
    if (strlen(chains) == 0) return NULL;

    freesasa_structure *new_p = freesasa_structure_new();
    
    if (!new_p) {
        mem_fail();
        return NULL;
    }
    
    new_p->model = p->model;

    for (int i = 0; i < p->number_atoms; ++i) {
        struct atom *ai = p->a[i];
        char c = ai->chain_label;
        if (strchr(chains,c) != NULL) {
            const double *v = freesasa_coord_i(p->xyz,i);
            int res = freesasa_structure_add_atom(new_p, ai->atom_name,
                                                  ai->res_name, ai->res_number,
                                                  c, v[0], v[1], v[2]);
            if (res == FREESASA_FAIL) {
                fail_msg("");
                freesasa_structure_free(new_p);
                return NULL;
            }
        }
    }
    if (new_p->number_atoms == 0) {
        freesasa_structure_free(new_p);
        new_p = NULL;
    }
    return new_p;
}

const char *
freesasa_structure_chain_labels(const freesasa_structure *structure)
{
    assert(structure);
    return structure->chains;
}

const coord_t *
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

const char *
freesasa_structure_atom_name(const freesasa_structure *p,
                             int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i]->atom_name;
}

const char*
freesasa_structure_atom_res_name(const freesasa_structure *p,
                                 int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i]->res_name;
}

const char*
freesasa_structure_atom_res_number(const freesasa_structure *p,
                                   int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i]->res_number;
}

char
freesasa_structure_atom_chain(const freesasa_structure *p,
                              int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i]->chain_label;
}
const char*
freesasa_structure_atom_symbol(const freesasa_structure *p,
                               int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i]->symbol;
}

const char*
freesasa_structure_atom_descriptor(const freesasa_structure *p,
                                   int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i]->descriptor;
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
                   const freesasa_structure *p)
{
    assert(p);
    assert(output);
    assert(result);
    assert(result->sasa);

    const double* values = result->sasa;
    const double* radii = p->radius;
    char buf[PDB_LINE_STRL+1], buf2[6];
    int n = freesasa_structure_n(p);
    if (p->model > 0) fprintf(output,"MODEL     %4d\n",p->model);
    else fprintf(output,             "MODEL        1\n");
    // Write ATOM entries
    for (int i = 0; i < n; ++i) {
        if (p->a[i]->line == NULL) {
            return freesasa_fail("in %s(): PDB input not valid or not present.",
                                 __func__);
        }
        strncpy(buf,p->a[i]->line,PDB_LINE_STRL);
        sprintf(&buf[54],"%6.2f%6.2f",radii[i],values[i]);
        fprintf(output,"%s\n",buf);
    }
    // Write TER  and ENDMDL lines
    errno = 0;
    strncpy(buf2,&buf[6],5);
    buf2[5]='\0';
    fprintf(output,"TER   %5d     %4s %c%4s\nENDMDL\n",
            atoi(buf2)+1, p->a[n-1]->res_name,
            p->a[n-1]->chain_label, p->a[n-1]->res_number);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

int
freesasa_structure_model(const freesasa_structure *structure)
{
    return structure->model;
}

const double *
freesasa_structure_coord_array(const freesasa_structure *structure)
{
    return freesasa_coord_all(structure->xyz);
}

const double *
freesasa_structure_radius(const freesasa_structure *structure)
{
    assert(structure);
    return structure->radius;
}

void
freesasa_structure_set_radius(freesasa_structure *structure, const double* radii)
{
    assert(structure);
    assert(radii);
    memcpy(structure->radius, radii, structure->number_atoms*sizeof(double));
}
