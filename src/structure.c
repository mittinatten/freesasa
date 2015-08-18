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

extern int freesasa_fail(const char *format, ...);
extern int freesasa_warn(const char *format, ...);

typedef struct {
    char res_name[PDB_ATOM_RES_NAME_STRL+1];
    char res_number[PDB_ATOM_RES_NUMBER_STRL+1];
    char atom_name[PDB_ATOM_NAME_STRL+1];
    char descriptor[FREESASA_STRUCTURE_DESCRIPTOR_STRL];
    char line[PDB_LINE_STRL+1];
    char chain_label;
} atom_t;

struct freesasa_structure {
    atom_t *a;
    freesasa_coord *xyz;
    int number_atoms;
    int number_residues;
    int number_chains;
    int *res_first_atom;
    char **res_desc;
};

freesasa_structure* freesasa_structure_new()
{
    freesasa_structure *p = malloc(sizeof(freesasa_structure));
    assert(p);
    p->number_atoms = 0;
    p->number_residues = 0;
    p->number_chains = 0;
    p->a = malloc(sizeof(atom_t));
    assert(p->a);
    p->xyz = freesasa_coord_new();
    p->res_first_atom = malloc(sizeof(int));
    assert(p->res_first_atom);
    p->res_desc = malloc(sizeof(char*));
    assert(p->res_desc);
    return p;
}

void freesasa_structure_free(freesasa_structure *p)
{
    if (p == NULL) return;
    if (p->a) free(p->a);
    if (p->xyz) freesasa_coord_free(p->xyz);
    if (p->res_first_atom) free(p->res_first_atom);
    if (p->res_desc) {
        for (int i = 0; i < p->number_residues; ++i) {
            if (p->res_desc[i]) free(p->res_desc[i]);
        }
        free(p->res_desc);
    }
    free(p);
}

/* returns alt_label if there is any */
char freesasa_structure_get_pdb_atom(atom_t *a, double *xyz, const char *line)
{
    assert(a); assert(xyz); assert(line);
    a->chain_label = freesasa_pdb_get_chain_label(line);
    freesasa_pdb_get_coord(xyz, line);
    freesasa_pdb_get_atom_name(a->atom_name, line);
    freesasa_pdb_get_res_name(a->res_name, line);
    freesasa_pdb_get_res_number(a->res_number, line);
    return freesasa_pdb_get_alt_coord_label(line);
}

freesasa_structure* freesasa_structure_from_pdb(FILE *pdb_file,
                                                int include_hetatm)
{
    assert(pdb_file);
    freesasa_structure *p = freesasa_structure_new();
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    char the_alt = ' ';

    assert(p);

    while (getline(&line, &len, pdb_file) != -1) {
        if (strncmp("ATOM",line,4)==0 ||
            ( (include_hetatm == 1) &&
              (strncmp("HETATM",line,6) == 0) )
            ) {
            double v[3];
            atom_t a;
            char alt;
            if (freesasa_pdb_ishydrogen(line)) continue;
            alt = freesasa_structure_get_pdb_atom(&a,v,line);
            if ((alt != ' ' && the_alt == ' ') || (alt == ' ')) {
                the_alt = alt;
            } else if (alt != ' ' && alt != the_alt) {
                continue;
            }

            freesasa_structure_add_atom(p,a.atom_name,a.res_name,a.res_number,
                                        a.chain_label, v[0], v[1], v[2]);
            strncpy(p->a[p->number_atoms-1].line,line,PDB_LINE_STRL);
        }
        if (strncmp("ENDMDL",line,4)==0) {
            if (p->number_atoms == 0) {
                freesasa_fail("input had ENDMDL before first ATOM entry.");
                freesasa_structure_free(p);
                p = NULL;
            }
            break;
        }
    }
    free(line);
    if (p->number_atoms == 0) {
        freesasa_fail("input had no ATOM entries.");
        freesasa_structure_free(p);
        return NULL;
    }
    return p;
}

static void freesasa_structure_alloc_one(freesasa_structure *p)
{
    assert(p);
    int na = ++p->number_atoms;
    p->a = realloc(p->a,sizeof(atom_t)*na);
    assert(p->a);
}

int freesasa_structure_add_atom(freesasa_structure *p,
                                const char *atom_name,
                                const char *residue_name,
                                const char *residue_number,
                                char chain_label,
                                double x, double y, double z)
{
    assert(p);
    assert(atom_name); assert(residue_name); assert(residue_number);

    // check input for consistency
    int na;
    atom_t *a;
    int validity = freesasa_classify_validate_atom(residue_name,atom_name);

    if (validity != FREESASA_SUCCESS) {
        return freesasa_warn("Skipping atom '%s' in residue '%s'",
                             atom_name,residue_name);
    }

    // allocate memory, increase number of atoms counter
    freesasa_structure_alloc_one(p);
    na = p->number_atoms;
    assert(na > 0);
    freesasa_coord_append_xyz(p->xyz,&x,&y,&z,1);
    a = &p->a[na-1];
    strcpy(a->atom_name,atom_name);
    strcpy(a->res_name,residue_name);
    strcpy(a->res_number,residue_number);
    sprintf(a->descriptor,"%c %s %s %s",chain_label,residue_number,residue_name,atom_name);
    a->chain_label = chain_label;
    a->line[0] = '\0';

    /* here we assume atoms are ordered sequentially, i.e. chains are
       not mixed in input: if two sequential atoms have different
       residue numbers or chain labels a new residue or new chain is
       assumed to begin */
    if (p->number_chains == 0)
        ++p->number_chains;
    if (na > 1 && chain_label != p->a[na-2].chain_label)
        ++p->number_chains;

    if (p->number_residues == 0) {
        ++p->number_residues;
        p->res_first_atom[0] = 0;
        p->res_desc[0] = malloc(FREESASA_STRUCTURE_DESCRIPTOR_STRL);
        assert(p->res_desc[0]);
        sprintf(p->res_desc[0],"%c %s %s",
                chain_label,residue_number,residue_name);
    }
    if (na > 1 && strcmp(residue_number,p->a[na-2].res_number)) {
        int naa = ++p->number_residues;
        p->res_first_atom = realloc(p->res_first_atom, sizeof(int)*naa);
        assert(p->res_first_atom);
        p->res_first_atom[naa-1] = na-1;
        p->res_desc = realloc(p->res_desc, sizeof(char*)*naa);
        assert(p->res_desc);
        p->res_desc[naa-1] = malloc(FREESASA_STRUCTURE_DESCRIPTOR_STRL);
        assert(p->res_desc[naa-1]);
        sprintf(p->res_desc[naa-1],"%c %s %s",
                chain_label,residue_number,residue_name);
    }
    return FREESASA_SUCCESS;
}

void freesasa_structure_r(double *r,
                          const freesasa_structure *p,
                          double (*atom2radius)(const char *res_name,
                                                const char *atom_name))
{
    assert(r); assert(p);
    for (int i = 0; i < p->number_atoms; ++i) {
        r[i] = atom2radius(p->a[i].res_name, p->a[i].atom_name);
    }
}

void freesasa_structure_r_def(double *r, const freesasa_structure *p)
{
    freesasa_structure_r(r,p,freesasa_classify_radius);
}

const freesasa_coord* freesasa_structure_xyz(const freesasa_structure *p)
{
    assert(p);
    return p->xyz;
}

int freesasa_structure_n(const freesasa_structure *p)
{
    assert(p);
    return p->number_atoms;
}

int freesasa_structure_n_residues(const freesasa_structure *p)
{
    assert(p);
    return p->number_residues;
}

const char* freesasa_structure_atom_name(const freesasa_structure *p,
                                         int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].atom_name;
}

const char* freesasa_structure_atom_res_name(const freesasa_structure *p,
                                             int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].res_name;
}

const char* freesasa_structure_atom_res_number(const freesasa_structure *p,
                                               int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].res_number;
}

char freesasa_structure_atom_chain(const freesasa_structure *p,
                                   int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].chain_label;
}

const char* freesasa_structure_atom_descriptor(const freesasa_structure *p, int i)
{
    assert(p);
    assert(i < p->number_atoms && i >= 0);
    return p->a[i].descriptor;
}

int freesasa_structure_residue_atoms(const freesasa_structure *s, int r_i, int *first, int *last)
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

const char* freesasa_structure_residue_descriptor(const freesasa_structure *s, int r_i)
{
    assert(s);
    assert(r_i < s->number_residues);
    return s->res_desc[r_i];
}

int freesasa_structure_write_pdb_bfactors(const freesasa_structure *p,
                                          FILE *output,
                                          const double *values,
                                          const double *radius)
{
    assert(p);
    assert(output);
    assert(values);
    assert(radius);

    // Write ATOM entries
    char buf[PDB_LINE_STRL+1], buf2[6];
    int n = freesasa_structure_n(p);
    for (int i = 0; i < n; ++i) {
        if (p->a[i].line[0] == '\0') {
            return freesasa_fail("%s: PDB input not valid or not present.",
                                 __func__);
        }
        strncpy(buf,p->a[i].line,PDB_LINE_STRL);
        sprintf(&buf[54],"%6.2f%6.2f",radius[i],values[i]);
        errno = 0;
        if (fprintf(output,"%s\n",buf) < 0)
            return freesasa_fail("%s: %s", __func__, strerror(errno));
    }
    // Write TER line
    errno = 0;
    strncpy(buf2,&buf[6],5);
    buf2[5]='\0';
    if (fprintf(output,"TER   %5d     %4s %c%4s\n",
                atoi(buf2)+1, p->a[n-1].res_name,
                p->a[n-1].chain_label, p->a[n-1].res_number) < 0) {
        freesasa_fail("%s: %s", __func__, strerror(errno));
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;
}
