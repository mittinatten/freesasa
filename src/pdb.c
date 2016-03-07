#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "freesasa_internal.h"
#include "pdb.h"

static inline int
pdb_line_check(const char *line,int len)
{
    assert(line);
    if (strlen(line) < 6) return FREESASA_FAIL;
    if (! strncmp(line,"ATOM",4) &&
        ! strncmp(line,"HETATM",6)) {
        return FREESASA_FAIL;
    }
    if (strlen(line) < len) {
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;
}

int
freesasa_pdb_get_models(FILE* pdb,
                        struct file_range** ranges)
{
    assert(pdb != NULL);
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    int n = 0, n_end = 0, error = 0;
    long last_pos = ftell(pdb);
    struct file_range *it = NULL, *itb;

    while (getline(&line, &len, pdb) != -1) {
        if (strncmp("MODEL",line,5)==0) {
            ++n;
            itb = it;
            it = realloc(it, sizeof(struct file_range)*n);
            if (!it) {
                free(itb);
                error = mem_fail();
                break;
            }
            it[n-1].begin = last_pos;
        }
        if (strncmp("ENDMDL",line,6)==0) {
            ++n_end;
            if (n != n_end) {
                error = freesasa_fail("in %s(): Mismatch between MODEL and ENDMDL "
                                      "in input\n",__func__);
                break;
            }
            it[n-1].end = ftell(pdb);
        }
        last_pos = ftell(pdb);
    }
    free(line);
    if (n == 0) { // when there are no models, the whole file is the model
        free(it);
        it = NULL;
    }
    if (error == FREESASA_FAIL) {
        free(it);
        *ranges = NULL;
        return FREESASA_FAIL;
    }
    *ranges = it;
    return n;
}

int
freesasa_pdb_get_chains(FILE *pdb,
                        struct file_range model,
                        struct file_range **ranges,
                        int options)
{
    assert(pdb);
    assert(ranges);
    // it is assumed that 'model' is valid for 'pdb'

    int n_chains = 0;
    size_t len = PDB_LINE_STRL;
    char *line = NULL;
    struct file_range *chains = NULL, *chb;
    char last_chain = '\0';
    long last_pos = model.begin;
    *ranges = NULL;

    // for each model, find file ranges for each chain, store them
    // in the dynamically growing array chains
    fseek(pdb,model.begin,SEEK_SET);
    while (getline(&line, &len, pdb) != -1 &&
           ftell(pdb) < model.end ) {
        if (strncmp("ATOM",line,4)==0 || ( (options & FREESASA_INCLUDE_HETATM) &&
                                           (strncmp("HETATM",line,6) == 0) ) ) {
            char chain = freesasa_pdb_get_chain_label(line);
            if (chain != last_chain) {
                if (n_chains > 0) chains[n_chains-1].end = last_pos;
                ++n_chains;
                chb = chains;
                chains = realloc(chains,sizeof(struct file_range)*n_chains);
                if (!chains) {
                    free(chb);
                    free(line);
                    return mem_fail();
                }
                chains[n_chains-1].begin = last_pos;
                last_chain = chain;
            }
        }
        last_pos = ftell(pdb);
    }
    free(line);

    if (n_chains > 0) {
        chains[n_chains-1].end = last_pos;
        chains[0].begin = model.begin; //preserve model info
        *ranges = chains;
    } else {
        *ranges = NULL;
    }
    return n_chains;
}


int
freesasa_pdb_get_atom_name(char *name,
                           const char *line)
{
    assert(name);
    assert(line);
    if (pdb_line_check(line,PDB_ATOM_NAME_STRL+12) == FREESASA_FAIL) {
        name[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(name,line+12,PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
    return FREESASA_SUCCESS;
}

int
freesasa_pdb_get_res_name(char *name,
                          const char *line)
{
    assert(name);
    assert(line);
    if (pdb_line_check(line,PDB_ATOM_RES_NAME_STRL+17) == FREESASA_FAIL) {
        name[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(name, line+17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
    return FREESASA_SUCCESS;
}

int
freesasa_pdb_get_coord(double *xyz,
                       const char *line)
{
    assert(xyz);
    assert(line);
    if (pdb_line_check(line,54) == FREESASA_FAIL) {
        return FREESASA_FAIL;
    }
    if (sscanf(line+30, "%lf%lf%lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
        return freesasa_fail("Could not read coordinates from line '%s'",line);
    }
    return FREESASA_SUCCESS;
}

int
freesasa_pdb_get_res_number(char *number,
                            const char* line)
{
    assert(number);
    assert(line);
    if (pdb_line_check(line,PDB_ATOM_RES_NUMBER_STRL+22) == FREESASA_FAIL) {
        number[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(number, line+22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
    return FREESASA_SUCCESS;
}
char
freesasa_pdb_get_chain_label(const char* line)
{
    assert(line);
    if (pdb_line_check(line,21) == FREESASA_FAIL) return '\0';
    return line[21];
}

char 
freesasa_pdb_get_alt_coord_label(const char* line)
{
    assert(line);
    if (pdb_line_check(line,16) == FREESASA_FAIL) return '\0';
    return line[16];
}

int
freesasa_pdb_get_symbol(char *symbol,
                        const char* line)
{
    assert(line);
    if (pdb_line_check(line,76+PDB_ATOM_SYMBOL_STRL) == FREESASA_FAIL) {
        symbol[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(symbol,line+76,2);
    symbol[2] = '\0';
    return FREESASA_SUCCESS;
}

int
freesasa_pdb_ishydrogen(const char* line)
{
    assert(line);
    if (pdb_line_check(line,13) == FREESASA_FAIL) return FREESASA_FAIL;
    //hydrogen
    if (line[12] == 'H' || line[13] == 'H') return 1;
    //hydrogen
    if (line[12] == 'D' || line[13] == 'D') return 1;
    return 0;
}
