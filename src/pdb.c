#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <errno.h>
#include <stdlib.h>

#include "freesasa_internal.h"
#include "pdb.h"

/* len >= 6 */
static inline int
pdb_line_check(const char *line, size_t len)
{
    assert(line);
    if (len < 6) return FREESASA_FAIL;
    if (strlen(line) < len) return FREESASA_FAIL;
    if (strncmp("ATOM", line, 4) != 0 &&
        strncmp("HETATM", line, 6) != 0) {
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;
}

/**
    Extracts a double from the line of maximum width characters, to
    allow checking for empty fields (instead of just reading the first
    float that comes along.
 */
static inline int
pdb_get_double(const char *line, size_t width, double *val)
{
    /* allow truncated lines */
    char buf[PDB_LINE_STRL];
    float tmp;

    if (strlen(line) < width) width = strlen(line);

    memcpy(buf, line, width);
    buf[width] = '\0';

    if (sscanf(buf, "%f", &tmp) == 1) {
        *val = tmp;
        return FREESASA_SUCCESS;
    }

    return FREESASA_FAIL;
}

int freesasa_pdb_get_models(FILE *pdb,
                            struct file_range **ranges)
{
    char line[PDB_MAX_LINE_STRL];
    int n = 0, n_end = 0, error = 0;
    long last_pos = ftell(pdb);
    struct file_range *it = NULL, *itb;

    assert(pdb != NULL);

    while (fgets(line, PDB_MAX_LINE_STRL, pdb) != NULL) {
        if (strncmp("MODEL", line, 5) == 0) {
            ++n;
            itb = it;
            it = realloc(it, sizeof(struct file_range) * n);
            if (!it) {
                free(itb);
                error = mem_fail();
                break;
            }
            it[n - 1].begin = last_pos;
        }
        if (strncmp("ENDMDL", line, 6) == 0) {
            ++n_end;
            if (n != n_end) {
                error = fail_msg("mismatch between MODEL and ENDMDL in input");
                break;
            }
            it[n - 1].end = ftell(pdb);
        }
        last_pos = ftell(pdb);
    }
    if (n == 0) { /* when there are no models, the whole file is the model */
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

int freesasa_pdb_get_chains(FILE *pdb,
                            struct file_range model,
                            struct file_range **ranges,
                            int options)
{
    /* it is assumed that 'model' is valid for 'pdb' */

    int n_chains = 0;
    char line[PDB_MAX_LINE_STRL];
    struct file_range *chains = NULL, *chb;
    char last_chain = '\0';
    long last_pos = model.begin;

    assert(pdb);
    assert(ranges);

    *ranges = NULL;

    /* for each model, find file ranges for each chain, store them
       in the dynamically growing array chains */
    fseek(pdb, model.begin, SEEK_SET);
    while (fgets(line, PDB_MAX_LINE_STRL, pdb) != NULL &&
           ftell(pdb) < model.end) {
        if (strncmp("ATOM", line, 4) == 0 || ((options & FREESASA_INCLUDE_HETATM) &&
                                              (strncmp("HETATM", line, 6) == 0))) {
            char chain = freesasa_pdb_get_chain_label(line);
            if (chain != last_chain) {
                if (n_chains > 0) chains[n_chains - 1].end = last_pos;
                ++n_chains;
                chb = chains;
                chains = realloc(chains, sizeof(struct file_range) * n_chains);
                if (!chains) {
                    free(chb);
                    return mem_fail();
                }
                chains[n_chains - 1].begin = last_pos;
                last_chain = chain;
            }
        }
        last_pos = ftell(pdb);
    }

    if (n_chains > 0) {
        chains[n_chains - 1].end = last_pos;
        chains[0].begin = model.begin; /* preserve model info */
        *ranges = chains;
    } else {
        *ranges = NULL;
    }
    return n_chains;
}

int freesasa_pdb_get_atom_name(char *name,
                               const char *line)
{
    assert(name);
    assert(line);
    if (pdb_line_check(line, PDB_ATOM_NAME_STRL + 12) == FREESASA_FAIL) {
        name[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(name, line + 12, PDB_ATOM_NAME_STRL);
    name[PDB_ATOM_NAME_STRL] = '\0';
    return FREESASA_SUCCESS;
}

int freesasa_pdb_get_res_name(char *name,
                              const char *line)
{
    assert(name);
    assert(line);
    if (pdb_line_check(line, PDB_ATOM_RES_NAME_STRL + 17) == FREESASA_FAIL) {
        name[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(name, line + 17, PDB_ATOM_RES_NAME_STRL);
    name[PDB_ATOM_RES_NAME_STRL] = '\0';
    return FREESASA_SUCCESS;
}

int freesasa_pdb_get_coord(double *xyz,
                           const char *line)
{
    int n_coord = 24; /* 54-30+1 */
    char coord_section[25];

    assert(xyz);
    assert(line);

    if (pdb_line_check(line, 54) == FREESASA_FAIL) {
        return FREESASA_FAIL;
    }

    strncpy(coord_section, line + 30, n_coord);
    coord_section[n_coord] = '\0';

    if (sscanf(coord_section, "%lf%lf%lf", &xyz[0], &xyz[1], &xyz[2]) != 3) {
        return fail_msg("could not read coordinates from line '%s'", line);
    }

    return FREESASA_SUCCESS;
}

int freesasa_pdb_get_res_number(char *number,
                                const char *line)
{
    assert(number);
    assert(line);
    if (pdb_line_check(line, PDB_ATOM_RES_NUMBER_STRL + 22) == FREESASA_FAIL) {
        number[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(number, line + 22, PDB_ATOM_RES_NUMBER_STRL);
    number[PDB_ATOM_RES_NUMBER_STRL] = '\0';
    return FREESASA_SUCCESS;
}
char freesasa_pdb_get_chain_label(const char *line)
{
    assert(line);
    if (pdb_line_check(line, 21) == FREESASA_FAIL) return '\0';
    return line[21];
}

char freesasa_pdb_get_alt_coord_label(const char *line)
{
    assert(line);
    if (pdb_line_check(line, 16) == FREESASA_FAIL) return '\0';
    return line[16];
}

int freesasa_pdb_get_symbol(char *symbol,
                            const char *line)
{
    assert(line);
    if (pdb_line_check(line, 76 + PDB_ATOM_SYMBOL_STRL) == FREESASA_FAIL) {
        symbol[0] = '\0';
        return FREESASA_FAIL;
    }
    strncpy(symbol, line + 76, 2);
    symbol[2] = '\0';
    return FREESASA_SUCCESS;
}

int freesasa_pdb_get_occupancy(double *occ,
                               const char *line)
{
    assert(line);
    /* allow truncated lines */
    if (pdb_line_check(line, 55) == FREESASA_SUCCESS)
        return pdb_get_double(line + 54, 6, occ);
    return FREESASA_FAIL;
}

int freesasa_pdb_get_bfactor(double *bfac,
                             const char *line)
{
    assert(line);
    /* allow truncated lines */
    if (pdb_line_check(line, 61) == FREESASA_SUCCESS)
        return pdb_get_double(line + 60, 6, bfac);
    return FREESASA_FAIL;
}

int freesasa_pdb_ishydrogen(const char *line)
{
    assert(line);
    char symbol[PDB_ATOM_SYMBOL_STRL + 1];
    freesasa_pdb_get_symbol(symbol, line);

    if (pdb_line_check(line, 13) == FREESASA_FAIL) return FREESASA_FAIL;

    /* Check symbol first */
    if (strncmp(symbol, " H", 2) == 0) return 1;
    if (strncmp(symbol, " D", 2) == 0) return 1;
    /* If the symbol is not blank and not H or D */
    if (!(strncmp(symbol, "  ", 2) == 0)) return 0;

    /* When symbol is missing */
    /* Cover elements such as Cd, Nd, Th, etc. ("CD  ", "ND  ") */
    if (!(line[12] == ' ' || (line[12] >= '1' && line[12] <= '9'))) return 0;
    /* Hydrogen, atom name "H**" or " H**" */
    if (line[12] == 'H' || line[13] == 'H') return 1;
    /* Deuterium */
    if (line[12] == 'D' || line[13] == 'D') return 1;
    return 0;
}

static int
write_pdb_impl(FILE *output,
               freesasa_node *structure)
{
    char buf[PDB_LINE_STRL + 1], buf2[6];
    int model;
    double radius;
    const char *line = NULL;
    freesasa_node *chain = NULL, *residue = NULL, *atom = NULL;
    const freesasa_nodearea *area = NULL;
    const char *last_res_name = NULL, *last_res_number = NULL, *last_chain = NULL;

    assert(freesasa_node_type(structure) == FREESASA_NODE_STRUCTURE);

    model = freesasa_node_structure_model(structure);
    if (model > 0)
        fprintf(output, "MODEL     %4d\n", model);
    else
        fprintf(output, "MODEL        1\n");

    chain = freesasa_node_children(structure);

    /* Write ATOM entries */
    while (chain) {
        residue = freesasa_node_children(chain);
        while (residue) {
            atom = freesasa_node_children(residue);
            while (atom) {
                line = freesasa_node_atom_pdb_line(atom);
                area = freesasa_node_area(atom);
                radius = freesasa_node_atom_radius(atom);

                if (line == NULL) {
                    return fail_msg("PDB input not valid or not present");
                }

                strncpy(buf, line, PDB_LINE_STRL);
                sprintf(&buf[54], "%6.2f%6.2f", radius, area->total);
                fprintf(output, "%s\n", buf);

                atom = freesasa_node_next(atom);
            }
            last_res_name = freesasa_node_name(residue);
            last_res_number = freesasa_node_residue_number(residue);
            residue = freesasa_node_next(residue);
        }
        last_chain = freesasa_node_name(chain);
        chain = freesasa_node_next(chain);
    }

    /* Write TER and ENDMDL lines */
    strncpy(buf2, &buf[6], 5);
    buf2[5] = '\0';
    fprintf(output, "TER   %5d     %4s %c%5s\nENDMDL\n",
            atoi(buf2) + 1, last_res_name, last_chain[0], last_res_number);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }

    return FREESASA_SUCCESS;
}

int freesasa_write_pdb(FILE *output,
                       freesasa_node *root)
{
    freesasa_node *result = freesasa_node_children(root), *structure;

    assert(output);
    assert(root);
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    fprintf(output, "REMARK 999 This PDB file was generated by %s.\n", freesasa_string);
    fprintf(output, "REMARK 999 In the ATOM records temperature factors have been\n"
                    "REMARK 999 replaced by the SASA of the atom, and the occupancy\n"
                    "REMARK 999 by the radius used in the calculation.\n");

    while (result) {
        structure = freesasa_node_children(result);
        while (structure) {
            if (write_pdb_impl(output, structure) == FREESASA_FAIL) {
                return fail_msg("");
            }
            structure = freesasa_node_next(structure);
        }
        result = freesasa_node_next(result);
    }

    return FREESASA_SUCCESS;
}

#if USE_CHECK
#include <check.h>
#include <math.h>

START_TEST(test_pdb)
{
    double v;
    ck_assert(pdb_line_check("", 6) == FREESASA_FAIL);
    ck_assert(pdb_line_check("ATOM", 4) == FREESASA_FAIL);
    ck_assert(pdb_line_check("ATOM", 6) == FREESASA_FAIL);
    ck_assert(pdb_line_check("HETAT", 6) == FREESASA_FAIL);
    ck_assert(pdb_line_check("HETAT ", 6) == FREESASA_FAIL);
    ck_assert(pdb_line_check("BLA BLA BLA", 10) == FREESASA_FAIL);
    ck_assert(pdb_line_check("BLA BLA BLA", 11) == FREESASA_FAIL);
    ck_assert(pdb_line_check("BLA BLA BLA", 12) == FREESASA_FAIL);

    /* these will pass, although they would be useless */
    ck_assert(pdb_line_check("ATOM  ", 6) == FREESASA_SUCCESS);
    ck_assert(pdb_line_check("HETATM", 6) == FREESASA_SUCCESS);

    /* a more likely type of error */
    ck_assert(pdb_line_check("HETATM", 7) == FREESASA_FAIL);

    /* The normal case */
    ck_assert(pdb_line_check("ATOM      1  N   MET A   1      27.340  "
                             "24.430   2.614  1.00  9.67           N  ",
                             80) == FREESASA_SUCCESS);

    v = 0;
    ck_assert(pdb_get_double("1.23", 4, &v) == FREESASA_SUCCESS);
    ck_assert(fabs(1.23 - v) < 1e-5);
    v = 0;
    ck_assert(pdb_get_double(" 1.23", 5, &v) == FREESASA_SUCCESS);
    ck_assert(fabs(1.23 - v) < 1e-5);
    v = 0;
    ck_assert(pdb_get_double("1.23", 10, &v) == FREESASA_SUCCESS);
    ck_assert(fabs(1.23 - v) < 1e-5);
    v = 0;
    ck_assert(pdb_get_double("1.23 4.56", 10, &v) == FREESASA_SUCCESS);
    ck_assert(fabs(1.23 - v) < 1e-5);

    ck_assert(pdb_get_double("  ", 10, &v) == FREESASA_FAIL);
    ck_assert(pdb_get_double("    1.23", 4, &v) == FREESASA_FAIL);
    ck_assert(pdb_get_double("abc", 10, &v) == FREESASA_FAIL);
    ck_assert(pdb_get_double("a 1.23", 6, &v) == FREESASA_FAIL);
}
END_TEST

TCase *
test_pdb_static()
{
    TCase *tc = tcase_create("pdb.c static");
    tcase_add_test(tc, test_pdb);

    return tc;
}

#endif /* USE_CHECK */
