#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <ctype.h>
#include <stdlib.h>

#include "freesasa_internal.h"
#include "pdb.h"
#include "selection.h"

#include "parser.h"

#include "lexer.h"

struct freesasa_selection {
    char *name;
    char *command;
    double area;
    int n_atoms;
};

struct selection {
    const char *name;
    int *atom;
    int size;
};

static const char *
e_str(expression_type e)
{
    switch (e) {
    case E_SELECTION:
        return "<selection>";
    case E_SYMBOL:
        return "symbol";
    case E_NAME:
        return "name";
    case E_RESN:
        return "resn";
    case E_RESI:
        return "resi";
    case E_CHAIN:
        return "chain";
    case E_ID:
        return "<id>";
    case E_NUMBER:
        return "<number>";
    case E_NEGNUM:
        return "<neg_number>";
    case E_AND:
        return "and";
    case E_OR:
        return "or";
    case E_NOT:
        return "not";
    case E_PLUS:
        return "< + >";
    case E_RANGE:
        return "< - >";
    case E_RANGE_OPEN_R:
        return "< - R >";
    case E_RANGE_OPEN_L:
        return "< - L >";
    }
    return NULL;
}

extern int freesasa_yyparse(expression **expression, yyscan_t scanner);

static expression *
expression_new()
{
    struct expression *e = malloc(sizeof(expression));

    if (e == NULL)
        mem_fail();
    else {
        e->type = E_SELECTION;
        e->left = NULL;
        e->right = NULL;
        e->value = NULL;
    }
    return e;
}

static void
expression_free(expression *expression)
{
    if (expression) {
        expression_free(expression->left);
        expression_free(expression->right);
        free(expression->value);
        free(expression);
    }
}

expression *
freesasa_selection_atom(expression_type type,
                        const char *val)
{
    expression *e = expression_new();
    int i, n;
    char *buf;

    assert(val);

    if (e != NULL) {
        if (type == E_NEGNUM) {
            n = strlen(val) + 2;
            buf = malloc(n);
            if (buf == NULL) {
                mem_fail();
                expression_free(e);
                return NULL;
            }
            sprintf(buf, "-%s", val);
            e->type = E_NUMBER;
            e->value = strdup(buf);
            free(buf);
        } else {
            e->type = type;
            e->value = strdup(val);
        }

        if (e->value == NULL) {
            mem_fail();
            expression_free(e);
            return NULL;
        }

        for (i = 0; i < strlen(e->value); ++i)
            e->value[i] = toupper(e->value[i]);
    }
    return e;
}

expression *
freesasa_selection_create(expression *selection,
                          const char *id)
{
    expression *e = expression_new();

    assert(id);

    if (e == NULL)
        expression_free(selection);
    else {
        e->type = E_SELECTION;
        e->left = selection;
        e->value = strdup(id);

        if (e->value == NULL) {
            mem_fail();
            expression_free(e);
            e = NULL;
        }
    }

    return e;
}

expression *
freesasa_selection_selector(expression_type type,
                            expression *list)
{
    expression *e = expression_new();

    if (e == NULL) {
        expression_free(list);
    } else {
        e->type = type;
        e->left = list;
    }

    return e;
}

expression *
freesasa_selection_operation(expression_type type,
                             expression *left,
                             expression *right)
{
    expression *e = expression_new();

    if (e == NULL) {
        expression_free(left);
        expression_free(right);
    } else {
        e->type = type;
        e->left = left;
        e->right = right;
    }

    return e;
}

/* for debugging */
static void
print_expr(FILE *output, const expression *e, int level)
{
    int i;
    fprintf(output, "\n");
    for (i = 0; i < level; ++i)
        fprintf(output, "  ");

    if (e == NULL) {
        fprintf(output, "()");
    } else {
        fprintf(output, "(%s ", e_str(e->type));
        if (e->value) fprintf(output, ": %s ", e->value);
        print_expr(output, e->left, level + 1);
        print_expr(output, e->right, level + 1);
        fprintf(output, ")");
    }
    fflush(output);
}

static expression *
get_expression(const char *selector)
{
    freesasa_yyscan_t scanner;
    YY_BUFFER_STATE state;
    int err;
    expression *expression = NULL;

    if (freesasa_yylex_init(&scanner)) {
        fail_msg("lexer failed");
    } else {
        state = freesasa_yy_scan_string(selector, scanner);
        err = freesasa_yyparse(&expression, scanner);
        if (err) {
            if (err == 1) fail_msg("parser failed");
            if (err == 2) mem_fail();
            expression_free(expression);
            expression = NULL;
        }
        freesasa_yy_delete_buffer(state, scanner);
        freesasa_yylex_destroy(scanner);
    }

    return expression;
}

static struct selection *
selection_new(int n)
{
    struct selection *selection = malloc(sizeof(struct selection));
    int i;

    if (selection == NULL) {
        mem_fail();
    } else {
        selection->size = n;
        selection->atom = malloc(sizeof(int) * n);

        if (selection->atom == NULL) {
            free(selection);
            mem_fail();
            selection = NULL;
        } else {
            for (i = 0; i < n; ++i)
                selection->atom[i] = 0;
        }
    }

    return selection;
}

static void
selection_free(struct selection *selection)
{
    if (selection) {
        free(selection->atom);
        free(selection);
    }
}

/* Looks for exact match between the atom-name and expr->value */
static int
match_name(const freesasa_structure *structure,
           const char *id,
           int i)
{
    char atom[PDB_ATOM_NAME_STRL + 1];
    sscanf(freesasa_structure_atom_name(structure, i), "%s", atom);
    if (strcmp(atom, id) == 0)
        return 1;
    return 0;
}

/** Looks for match of strlen(expr->value) first characters of atom-name and expr->value */
static int
match_symbol(const freesasa_structure *structure,
             const char *id,
             int i)
{
    char symbol[PDB_ATOM_SYMBOL_STRL + 1];
    sscanf(freesasa_structure_atom_symbol(structure, i), "%s", symbol);
    if (strcmp(id, symbol) == 0) {
        return 1;
    }
    return 0;
}

static int
match_resn(const freesasa_structure *structure,
           const char *id,
           int i)
{
    char resn[PDB_ATOM_RES_NAME_STRL + 1];
    sscanf(freesasa_structure_atom_res_name(structure, i), "%s", resn);
    if (strcmp(resn, id) == 0)
        return 1;
    return 0;
}

static int
match_resi(const freesasa_structure *structure,
           const char *id,
           int i)
{
    char resi[PDB_ATOM_RES_NUMBER_STRL + 1];
    sscanf(freesasa_structure_atom_res_number(structure, i), "%s", resi);
    if (strcmp(resi, id) == 0)
        return 1;
    return 0;
}

static int
match_chain(const freesasa_structure *structure,
            const char *id,
            int i)
{
    return id[0] == freesasa_structure_atom_chain(structure, i);
}

static void
select_id(expression_type parent_type,
          struct selection *selection,
          const freesasa_structure *structure,
          const char *id)
{
    int count = 0, match, i;

    assert(id);

    for (i = 0; i < selection->size; ++i) {
        match = 0;
        switch (parent_type) {
        case E_NAME:
            match = match_name(structure, id, i);
            break;
        case E_SYMBOL:
            match = match_symbol(structure, id, i);
            break;
        case E_RESN:
            match = match_resn(structure, id, i);
            break;
        case E_RESI:
            match = match_resi(structure, id, i);
            break;
        case E_CHAIN:
            match = match_chain(structure, id, i);
            break;
        default:
            assert(0);
            break;
        }
        if (match) selection->atom[i] = 1;
        count += match;
    }
    if (count == 0) freesasa_warn("Found no matches to %s '%s', typo?",
                                  e_str(parent_type), id);
}

static int
is_valid_id(int parent_type,
            const expression *expr)
{
    int type, n, warn, i;
    const char *val;

    assert(expr->type == E_NUMBER || expr->type == E_ID);

    type = expr->type;
    val = expr->value;

    switch (parent_type) {
    case E_NAME:
        if (strlen(val) > PDB_ATOM_NAME_STRL)
            return freesasa_warn("select: %s: atom name '%s' invalid (string too long), "
                                 "will be ignored",
                                 e_str(parent_type), val);
        break;
    case E_SYMBOL:
        if (type != E_ID)
            return freesasa_warn("select: %s: '%s' invalid (should be 1 or 2 letters, 'C', 'N', 'SE', etc), "
                                 "will be ignored",
                                 e_str(parent_type), val);
        if (strlen(val) > 2)
            return freesasa_warn("select: %s: '%s' invalid (element names have 1 or 2 characters), "
                                 "will be ignored",
                                 e_str(parent_type), val);
        break;
    case E_RESN:
        if (strlen(val) > PDB_ATOM_RES_NAME_STRL)
            return freesasa_warn("select: %s: '%s' invalid (string too long), "
                                 "will be ignored",
                                 e_str(parent_type), val);
        break;
    case E_RESI:
        if (type == E_ID) {
            /* these should have format 1, 2, 345, etc or 12A, 12B, etc. */
            n = strlen(val);
            if (n > PDB_ATOM_RES_NUMBER_STRL) {
                return freesasa_warn("select: %s: '%s' invalid (string too long), "
                                     "will be ignored",
                                     e_str(parent_type), val);
            } else {
                warn = 0;
                if (n == 1) ++warn;
                if (!warn && (toupper(val[n - 1]) < 'A' || toupper(val[n - 1]) > 'Z')) {
                    ++warn;
                }
                for (i = 0; !warn && i < n - 1; ++i) {
                    if (val[i] < '0' || val[i] > '9') {
                        ++warn;
                    }
                }
                if (warn) {
                    return freesasa_warn("select: %s: '%s' invalid, should either be "
                                         "number (1, 2, 3) or number with insertion code (1A, 1B, ...), "
                                         "will be ignored",
                                         e_str(parent_type), val);
                }
            }
        } else if (type != E_NUMBER) {
            return freesasa_warn("select: %s: '%s' invalid, will be ignored",
                                 e_str(parent_type), val);
        }
        break;
    case E_CHAIN:
        if (strlen(val) > 1)
            return freesasa_warn("select: %s: '%s' invalid (string too long), "
                                 "will be ignored",
                                 e_str(parent_type), val);
        break;
    default:
        assert(0);
        break;
    }
    return FREESASA_SUCCESS;
}

static int
select_range(expression_type range_type,
             expression_type parent_type,
             struct selection *selection,
             const freesasa_structure *structure,
             const expression *left,
             const expression *right)
{
    int lower, upper, i, j;

    assert(range_type == E_RANGE || range_type == E_RANGE_OPEN_L || range_type == E_RANGE_OPEN_R);
    assert(parent_type == E_RESI || parent_type == E_CHAIN);

    if (parent_type == E_RESI) { /* residues have integer numbering */
        if ((left && left->type != E_NUMBER) ||
            (right && right->type != E_NUMBER)) {
            return freesasa_warn("select: %s: range '%s-%s' invalid, needs to be two numbers, "
                                 "will be ignored",
                                 e_str(parent_type), left->value, right->value);
        }
    } else { /* chains can be numbered by both letters (common) and numbers (uncommon) */
        if (left->type != right->type ||
            (left->type == E_ID && (strlen(left->value) > 1 || strlen(right->value) > 1)))
            return freesasa_warn("select: %s: range '%s-%s' invalid, should be two letters (A-C) or numbers (1-5), "
                                 "will be ignored",
                                 e_str(parent_type), left->value, right->value);
    }
    if (range_type == E_RANGE_OPEN_L) {
        lower = atoi(freesasa_structure_atom_res_number(structure, 0));
        upper = atoi(right->value);
    } else if (range_type == E_RANGE_OPEN_R) {
        lower = atoi(left->value);
        upper = atoi(freesasa_structure_atom_res_number(structure, freesasa_structure_n(structure) - 1));
    } else if (left->type == E_NUMBER) {
        lower = atoi(left->value);
        upper = atoi(right->value);
    } else {
        lower = (int)left->value[0];
        upper = (int)right->value[0];
    }
    for (i = 0; i < selection->size; ++i) {
        if (parent_type == E_RESI)
            j = atoi(freesasa_structure_atom_res_number(structure, i));
        else
            j = (int)freesasa_structure_atom_chain(structure, i);
        if (j >= lower && j <= upper)
            selection->atom[i] = 1;
    }
    return FREESASA_SUCCESS;
}

static int
select_list(expression_type parent_type,
            struct selection *selection,
            const freesasa_structure *structure,
            const expression *expr)
{
    int resr, resl;
    expression *left, *right;

    if (expr == NULL)
        return fail_msg("NULL expression");

    left = expr->left;
    right = expr->right;

    switch (expr->type) {
    case E_PLUS:
        if (left == NULL || right == NULL)
            return fail_msg("NULL expression");
        resl = select_list(parent_type, selection, structure, left);
        resr = select_list(parent_type, selection, structure, right);
        if (resl == FREESASA_WARN || resr == FREESASA_WARN)
            return FREESASA_WARN;
        break;
    case E_RANGE:
        if (left == NULL || right == NULL)
            return fail_msg("NULL expression");
        return select_range(E_RANGE, parent_type, selection, structure, left, right);
    case E_RANGE_OPEN_L:
        if (left != NULL || right == NULL)
            return fail_msg("NULL expression");
        return select_range(E_RANGE_OPEN_L, parent_type, selection, structure, left, right);
    case E_RANGE_OPEN_R:
        if (left == NULL || right != NULL)
            return fail_msg("NULL expression");
        return select_range(E_RANGE_OPEN_R, parent_type, selection, structure, left, right);
    case E_ID:
    case E_NUMBER:
        if (is_valid_id(parent_type, expr) == FREESASA_SUCCESS)
            select_id(parent_type, selection, structure, expr->value);
        else
            return freesasa_warn("select: %s: '%s' invalid %s",
                                 e_str(parent_type), expr->value, e_str(expr->type));
        break;
    default:
        return freesasa_fail("select: parse error (expression: '%s %s')",
                             e_str(parent_type), e_str(expr->type));
    }
    return FREESASA_SUCCESS;
}

static int
selection_join(struct selection *target,
               const struct selection *s1,
               const struct selection *s2,
               int type)
{
    int n, i;

    if (s1 == NULL || s2 == NULL || target == NULL)
        return fail_msg("trying to join NULL selections");

    assert(s1->size == s2->size);
    assert(s1->size == target->size);

    n = target->size;

    switch (type) {
    case E_AND:
        for (i = 0; i < n; ++i)
            target->atom[i] = s1->atom[i] && s2->atom[i];
        break;
    case E_OR:
        for (i = 0; i < n; ++i)
            target->atom[i] = s1->atom[i] || s2->atom[i];
        break;
    default:
        assert(0);
    }

    return FREESASA_SUCCESS;
}

static int
selection_not(struct selection *s)
{
    int i;

    if (s == NULL) return fail_msg("NULL selection");

    for (i = 0; i < s->size; ++i) {
        s->atom[i] = !s->atom[i];
    }

    return FREESASA_SUCCESS;
}

/* Called recursively, the selection is built as we cover the expression tree */
static int
select_atoms(struct selection *selection,
             const expression *expr,
             const freesasa_structure *structure)
{
    int warn = 0, err = 0, n, ret;
    struct selection *sl, *sr;

    assert(selection);
    assert(structure);

    n = selection->size;

    /* this should only happen if memory allocation failed during parsing */
    if (expr == NULL) return fail_msg("NULL expression");

    switch (expr->type) {
    case E_SELECTION:
        assert(expr->value != NULL);
        selection->name = expr->value;
        return select_atoms(selection, expr->left, structure);
        break;
    case E_SYMBOL:
    case E_NAME:
    case E_RESN:
    case E_RESI:
    case E_CHAIN:
        return select_list(expr->type, selection, structure, expr->left);
        break;
    case E_AND:
    case E_OR: {
        sl = selection_new(n);
        sr = selection_new(n);

        if (sl != NULL && sr != NULL) {
            if ((ret = select_atoms(sl, expr->left, structure))) {
                if (ret == FREESASA_WARN) ++warn;
                if (ret == FREESASA_FAIL) ++err;
            }
            if ((ret = select_atoms(sr, expr->right, structure))) {
                if (ret == FREESASA_WARN) ++warn;
                if (ret == FREESASA_FAIL) ++err;
            }
            selection_join(selection, sl, sr, expr->type);
        } else {
            ++err;
        }

        selection_free(sl);
        selection_free(sr);

        if (err) return fail_msg("error joining selections");
        break;
    }
    case E_NOT: {
        ret = select_atoms(selection, expr->right, structure);
        if (ret == FREESASA_WARN) ++warn;
        if (ret == FREESASA_FAIL) return FREESASA_FAIL;
        if (selection_not(selection)) return FREESASA_FAIL;
        break;
    }
    case E_ID:
    case E_NUMBER:
    case E_PLUS:
    case E_RANGE:
        /* these four are handled by the RESN,SYMBOL,ETC */
    default:
        return fail_msg("parser error");
    }
    if (warn) return FREESASA_WARN;
    return FREESASA_SUCCESS;
}

static int
select_area_impl(const char *command,
                 char *name,
                 double *area,
                 const freesasa_structure *structure,
                 const freesasa_result *result)
{
    struct selection *selection = NULL;
    struct expression *expression = NULL;
    const int maxlen = FREESASA_MAX_SELECTION_NAME;
    double sasa = 0;
    int err = 0, warn = 0, n_atoms = 0, j, len;

    assert(name);
    assert(area);
    assert(command);
    assert(structure);
    assert(result);
    assert(freesasa_structure_n(structure) == result->n_atoms);

    *area = 0;
    name[0] = '\0';

    expression = get_expression(command);
    selection = selection_new(result->n_atoms);

    if (selection == NULL) {
        return fail_msg("");
    }

    if (expression != NULL && selection != NULL) {
        switch (select_atoms(selection, expression, structure)) {
        case FREESASA_FAIL:
            err = 1;
            break;
        case FREESASA_WARN:
            warn = 1; /* proceed with calculation, print warning later */
        case FREESASA_SUCCESS: {
            for (j = 0; j < selection->size; ++j) {
                ++n_atoms;
                sasa += selection->atom[j] * result->sasa[j];
            }

            *area = sasa;
            len = strlen(selection->name);
            strncpy(name, selection->name, maxlen);

            if (len > maxlen) {
                name[maxlen] = '\0';
            } else {
                name[len] = '\0';
            }

            break;
        }
        default:
            assert(0);
        }
    } else {
        err = 1;
    }
    selection_free(selection);
    expression_free(expression);

    if (err)
        return fail_msg("problems parsing expression '%s'", command);
    if (warn)
        return freesasa_warn("in %s(): There were warnings", __func__);
    return n_atoms;
}

freesasa_selection *
freesasa_selection_alloc(const char *name, const char *command)
{
    freesasa_selection *selection = malloc(sizeof(freesasa_selection));

    if (selection == NULL) {
        mem_fail();
        return NULL;
    }

    selection->name = NULL;
    selection->command = NULL;
    selection->area = 0;
    selection->n_atoms = 0;

    selection->name = strdup(name);
    if (selection->name == NULL) {
        mem_fail();
        goto cleanup;
    }

    selection->command = strdup(command);
    if (selection->command == NULL) {
        mem_fail();
        goto cleanup;
    }

    return selection;

cleanup:
    freesasa_selection_free(selection);
    return NULL;
}

void freesasa_selection_free(freesasa_selection *selection)
{
    if (selection != NULL) {
        free(selection->name);
        free(selection->command);
        free(selection);
    }
}

freesasa_selection *
freesasa_selection_clone(const freesasa_selection *src)
{
    freesasa_selection *cpy = freesasa_selection_alloc(src->name, src->command);

    if (cpy == NULL) {
        fail_msg("");
        goto cleanup;
    }

    cpy->area = src->area;
    cpy->n_atoms = src->n_atoms;

    return cpy;

cleanup:
    freesasa_selection_free(cpy);
    return NULL;
}

const char *
freesasa_selection_name(const freesasa_selection *selection)
{
    assert(selection);
    return selection->name;
}

const char *
freesasa_selection_command(const freesasa_selection *selection)
{
    assert(selection);
    return selection->command;
}

double
freesasa_selection_area(const freesasa_selection *selection)
{
    assert(selection);
    return selection->area;
}

freesasa_selection *
freesasa_selection_new(const char *command,
                       const freesasa_structure *structure,
                       const freesasa_result *result)
{
    char name[FREESASA_MAX_SELECTION_NAME];
    double area;
    freesasa_selection *selection;
    int n_atoms;

    n_atoms = select_area_impl(command, name, &area, structure, result);

    if (n_atoms == FREESASA_FAIL) {
        fail_msg("");
        return NULL;
    }

    selection = freesasa_selection_alloc(name, command);
    if (selection == NULL) {
        mem_fail();
        return NULL;
    }

    selection->area = area;
    selection->n_atoms = n_atoms;

    return selection;
}

int freesasa_select_area(const char *command,
                         char *name,
                         double *area,
                         const freesasa_structure *structure,
                         const freesasa_result *result)
{
    int ret = select_area_impl(command, name, area, structure, result);
    if (ret >= 0) return FREESASA_SUCCESS;
    return ret;
}

int freesasa_selection_parse_error(expression *e,
                                   yyscan_t scanner,
                                   const char *msg)
{
    if (freesasa_get_verbosity() == FREESASA_V_DEBUG) print_expr(stderr, e, 0);
    if (freesasa_get_verbosity() == FREESASA_V_NORMAL) fprintf(stderr, "\n");
    return freesasa_fail(msg);
}

#if USE_CHECK
#include <check.h>

START_TEST(test_selection)
{
    struct selection *s1, *s2, *s3, *s4;
    static const expression empty_expression = {
        .right = NULL, .left = NULL, .value = NULL, .type = E_SELECTION};
    freesasa_structure *structure = freesasa_structure_new();
    expression r, l, e, e_symbol;

    freesasa_structure_add_atom(structure, " CA ", "ALA", "   1", 'A', 0, 0, 0);
    freesasa_structure_add_atom(structure, " O  ", "ALA", "   1", 'A', 10, 10, 10);

    s1 = selection_new(freesasa_structure_n(structure));
    s2 = selection_new(freesasa_structure_n(structure));
    s3 = selection_new(freesasa_structure_n(structure));
    s4 = selection_new(freesasa_structure_n(structure));

    r = l = e = e_symbol = empty_expression;
    e.type = E_PLUS;
    e.right = &r;
    e.left = &l;
    r.value = "C";
    r.type = E_ID;
    l.value = "O";
    l.type = E_ID;
    e_symbol.type = E_SYMBOL;
    e_symbol.left = &e;

    /* select_symbol */
    select_list(E_SYMBOL, s1, structure, &r);
    ck_assert_int_eq(s1->atom[0], 1);
    ck_assert_int_eq(s1->atom[1], 0);
    select_list(E_SYMBOL, s2, structure, &l);
    ck_assert_int_eq(s2->atom[0], 0);
    ck_assert_int_eq(s2->atom[1], 1);
    select_list(E_SYMBOL, s3, structure, &e);
    ck_assert_int_eq(s3->atom[0], 1);
    ck_assert_int_eq(s3->atom[1], 1);
    select_atoms(s4, &e_symbol, structure);
    ck_assert_int_eq(s4->atom[0], 1);
    ck_assert_int_eq(s4->atom[1], 1);

    /* selection_join */
    selection_join(s3, s1, s2, E_AND);
    ck_assert_int_eq(s3->atom[0], 0);
    ck_assert_int_eq(s3->atom[1], 0);
    selection_join(s3, s1, s2, E_OR);
    ck_assert_int_eq(s3->atom[0], 1);
    ck_assert_int_eq(s3->atom[1], 1);
    freesasa_set_verbosity(FREESASA_V_SILENT);
    ck_assert_int_eq(selection_join(NULL, s1, s2, E_OR), FREESASA_FAIL);
    ck_assert_int_eq(selection_join(s3, NULL, s1, E_OR), FREESASA_FAIL);
    ck_assert_int_eq(selection_join(NULL, NULL, NULL, E_OR), FREESASA_FAIL);

    /* selection_not */
    ck_assert_int_eq(selection_not(s3), FREESASA_SUCCESS);
    ck_assert_int_eq(s3->atom[0], 0);
    ck_assert_int_eq(s3->atom[1], 0);
    ck_assert_int_eq(selection_not(NULL), FREESASA_FAIL);
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

START_TEST(test_expression)
{
    int i;
    expression *e = get_expression("c1, symbol O+C");

    ck_assert_ptr_ne(e, NULL);
    ck_assert_int_eq(e->type, E_SELECTION);
    ck_assert_ptr_ne(e->left, NULL);
    ck_assert_ptr_eq(e->right, NULL);
    ck_assert_str_eq(e->value, "c1");
    ck_assert_int_eq(e->left->type, E_SYMBOL);
    ck_assert_ptr_ne(e->left->left, NULL);
    ck_assert_ptr_eq(e->left->right, NULL);
    ck_assert_int_eq(e->left->left->type, E_PLUS);
    ck_assert_int_eq(e->left->left->left->type, E_ID);
    ck_assert_int_eq(e->left->left->right->type, E_ID);
    ck_assert_str_eq(e->left->left->right->value, "C");
    ck_assert_str_eq(e->left->left->left->value, "O");
    for (i = E_SELECTION; i <= E_RANGE_OPEN_R; ++i) {
        ck_assert_ptr_ne(e_str(i), NULL);
    }
    ck_assert_ptr_eq(e_str(E_RANGE_OPEN_R + 1), NULL);
}
END_TEST

START_TEST(test_debug) /* this test just runs the debug output code to not get artificially low coverage */
{
    FILE *devnull = fopen("/dev/null", "w");
    expression *e = get_expression("c1, symbol O+C");
    print_expr(devnull, e, 0);
}
END_TEST

struct selection selection_dummy = {.size = 1, .name = NULL, .atom = NULL};

void *freesasa_selection_dummy_ptr = &selection_dummy;

int freesasa_wrap_select_atoms(struct selection *selection,
                               const expression *expr,
                               const freesasa_structure *structure)
{
    return select_atoms(selection, expr, structure);
}

TCase *
test_selection_static()
{
    TCase *tc = tcase_create("selection.c static");
    tcase_add_test(tc, test_selection);
    tcase_add_test(tc, test_expression);
    tcase_add_test(tc, test_debug);

    return tc;
}

#endif /* USE_CHECK */
