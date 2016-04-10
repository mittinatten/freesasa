#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include "selection.h"
#include "parser.h"
#include "lexer.h"
#include "freesasa.h"
#include "freesasa_internal.h"
#include "pdb.h"

struct selection {
    const char* name;
    int *atom;
    int size;
};


static const char*
e_str(expression_type e)
{
    switch(e) {
    case E_SELECTION: return "<selection>";
    case E_SYMBOL:    return "symbol";
    case E_NAME:      return "name";
    case E_RESN:      return "resn";
    case E_RESI:      return "resi";
    case E_CHAIN:     return "chain";
    case E_ID:        return "<id>";
    case E_NUMBER:    return "<number>";
    case E_AND:       return "and";
    case E_OR:        return "or";
    case E_NOT:       return "not";
    case E_PLUS:      return "< + >";
    case E_RANGE:     return "< - >";
    }
    return NULL;
}

extern int freesasa_yyparse(expression **expression, yyscan_t scanner);

static expression *
expression_new() 
{
    struct expression *e = malloc(sizeof(expression));

    if (e == NULL) {mem_fail(); return NULL; }
    
    e->type = E_SELECTION;
    e->left = NULL;
    e->right = NULL;
    e->value = NULL;
   
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
                        const char* val)
{
    assert(val);
    expression *e = expression_new();
    if (e == NULL) return NULL;

    e->type = type;
    e->value = strdup(val);

    if (e->value == NULL) {
        mem_fail();
        expression_free(e);
        return NULL;
    }

    for (int i = 0; i < strlen(val); ++i) e->value[i] = toupper(val[i]);
    
    return e;
}

expression *
freesasa_selection_create(expression *selection,
                          const char* id)
{
    assert(id);
    expression *e = expression_new();
    if (e == NULL) {
        expression_free(selection);
        return NULL;
    }

    e->type = E_SELECTION;
    e->left = selection;
    e->value = strdup(id);
    
    if (e->value == NULL) {
        mem_fail();
        expression_free(e);
        return NULL;
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
        return NULL;
    }
    e->type = type;
    e->left = list;
    
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
        return NULL;
    }
    e->type = type;
    e->left = left;
    e->right = right;

    return e;
}


static expression *
get_expression(const char *selector) 
{
   freesasa_yyscan_t scanner;
   YY_BUFFER_STATE state;
   int err;
   expression *expression = NULL;
   if (freesasa_yylex_init(&scanner)) {
       fail_msg("Lexer failed");
       return NULL;
   }
   state = freesasa_yy_scan_string(selector, scanner);
   err = freesasa_yyparse(&expression, scanner);
   if (err) {
       if (err == 1) fail_msg("Parser failed");
       if (err == 2) mem_fail();
       expression_free(expression);
   }
   freesasa_yy_delete_buffer(state, scanner);

   freesasa_yylex_destroy(scanner);

   if (err) return NULL;

   return expression;
}

static struct selection *
selection_new(int n)
{
    struct selection *selection = malloc(sizeof(struct selection));
    
    if (selection == NULL) { 
        mem_fail(); 
        return NULL; 
    }
    
    selection->size = n;
    selection->atom = calloc(n,sizeof(int));
    
    if (selection->atom == NULL) {
        free(selection);
        mem_fail();
        return NULL;
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

/** Looks for exact match between the atom-name and expr->value*/
static int
match_name(const freesasa_structure *structure,
           const char *id, 
           int i)
{
    char atom[PDB_ATOM_NAME_STRL+1];
    sscanf(freesasa_structure_atom_name(structure,i), "%s", atom);
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
    char symbol[PDB_ATOM_SYMBOL_STRL+1];
    sscanf(freesasa_structure_atom_symbol(structure,i), "%s", symbol);
    if (strcmp(id,symbol) == 0) {
        return 1;
    }
    return 0;
}

static int
match_resn(const freesasa_structure *structure,
           const char *id,
           int i)
{
    char resn[PDB_ATOM_RES_NAME_STRL+1];
    sscanf(freesasa_structure_atom_res_name(structure,i), "%s", resn);
    if (strcmp(resn, id) == 0)
        return 1;
    return 0;
}

static int
match_resi(const freesasa_structure *structure,
           const char *id, 
           int i)
{
    int resi = atoi(freesasa_structure_atom_res_number(structure,i));
    int e_resi = atoi(id);
    return resi == e_resi;
}

static int
match_chain(const freesasa_structure *structure,
            const char *id,
            int i)
{
    return id[0] == freesasa_structure_atom_chain(structure,i);
}

static void
select_id(expression_type parent_type,
          struct selection *selection,
          const freesasa_structure *structure,
          const char *id)
{
    assert(id);
    int count = 0;
    for (int i = 0; i < selection->size; ++i) {
        int match = 0;
        switch(parent_type) {
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
                                  e_str(parent_type),id);
}

static int
is_valid_id(int parent_type,
            const expression *expr)
{
    assert(expr->type == E_NUMBER || expr->type == E_ID);
    int type = expr->type;
    const char *val = expr->value;
    switch(parent_type) {
    case E_NAME:
        if (strlen(val) > PDB_ATOM_NAME_STRL)
            return freesasa_warn("select: Atom name '%s' invalid (string too long). "
                                 "Will be ignored",val);
        break;
    case E_SYMBOL:
        if (type != E_ID)
            return freesasa_warn("select: Symbol '%s' invalid (should be 1 or 2 letters, 'C', 'N', 'SE', etc). "
                                 "Will be ignored.",val);
        if (strlen(val) > 2)
            return freesasa_warn("select: Symbol '%s' invalid (element names have 1 or 2 characters). "
                                 "Will be ignored.",val);
        break;
    case E_RESN:
        if (strlen(val) > PDB_ATOM_RES_NAME_STRL)
            return freesasa_warn("select: Residue '%s' invalid (string too long). "
                                 "Will be ignored.",val);
        break;
    case E_RESI:
        if (type != E_NUMBER)
            return freesasa_warn("select: Residue number '%s' invalid (not a number). "
                                 "Will be ignored.",val);
        break;
    case E_CHAIN:
        if (strlen(val) > 1)
            return freesasa_warn("select: Chain label '%s' invalid (string too long). "
                                 "Will be ignored.",val);
        break;
    default:
        assert(0);
        break;
    }
    return FREESASA_SUCCESS;
}

static int
select_range(expression_type parent_type,
             struct selection *selection,
             const freesasa_structure *structure,
             const expression *left,
             const expression *right)
{
    assert(parent_type == E_RESI || parent_type == E_CHAIN);
    int lower, upper;
    if (parent_type == E_RESI) { // residues have integer numbering
        if (left->type != E_NUMBER || right->type != E_NUMBER) {
            return freesasa_warn("select: Range '%s-%s' invalid, needs to be two numbers (1-5). "
                                 "Will be ignored.",left->value,right->value);
        }
    } else { // chains can be numbered by both letters (common) and numbers (uncommon)
        if (left->type != right->type ||
            (left->type == E_ID && (strlen(left->value) > 1 || strlen(right->value) > 1)))
            return freesasa_warn("select: Chain range '%s-%s' invalid, needs to be two letters (A-C) or two numbers (1-5). "
                                 "Will be ignored.",left->value,right->value);
    }
    if (left->type == E_NUMBER) {
        lower = atoi(left->value);
        upper = atoi(right->value);
    } else {
        lower = (int)left->value[0];
        upper = (int)right->value[0];
    }
    for (int i = 0; i < selection->size; ++i) {
        int j;
        if (parent_type == E_RESI) j = atoi(freesasa_structure_atom_res_number(structure,i));
        else j = (int)freesasa_structure_atom_chain(structure,i);
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
    if (expr == NULL)
        return fail_msg("NULL expression.");
    int resr, resl;
    expression *left = expr->left, *right = expr->right;
    switch(expr->type) {
    case E_PLUS: 
        if (left == NULL || right == NULL) 
            return fail_msg("NULL expression.");
        resl = select_list(parent_type,selection,structure,left);
        resr = select_list(parent_type,selection,structure,right);
        if (resl == FREESASA_WARN || resr == FREESASA_WARN)
            return FREESASA_WARN;
        break;
    case E_RANGE:
        if (left == NULL || right == NULL) 
            return fail_msg("NULL expression.");
        return select_range(parent_type,selection,structure,left,right);
    case E_ID:
    case E_NUMBER:
        if (is_valid_id(parent_type,expr) == FREESASA_SUCCESS)
            select_id(parent_type,selection,structure,expr->value);
        else return FREESASA_WARN;
        break;
    default:
        return freesasa_fail("in %s(): parse error (expression: '%s %s')",
                             __func__,e_str(parent_type),e_str(expr->type));
    }
    return FREESASA_SUCCESS;
}

static int
selection_join(struct selection *target,
               const struct selection *s1,
               const struct selection *s2,
               int type)
{
    int n;
    if (s1 == NULL || s2 == NULL || target == NULL)
        return fail_msg("Trying to join NULL selections");
    
    assert(s1->size == s2->size);
    assert(s1->size == target->size);

    n = target->size;

    switch (type) {
    case E_AND:
        for (int i = 0; i < n; ++i)
            target->atom[i] = s1->atom[i] && s2->atom[i];
        break;
    case E_OR:
        for (int i = 0; i < n; ++i)
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
    if (s == NULL) return fail_msg("NULL selection");
    for (int i = 0; i < s->size; ++i) {
        s->atom[i] = ! s->atom[i];
    }
    return FREESASA_SUCCESS;
}

/* Called recursively, the selection is built as we cover the expression tree */
static int
select_atoms(struct selection* selection,
             const expression *expr,
             const freesasa_structure *structure)
{
    int warn = 0, err = 0, n = selection->size, ret;

    // this should only happen if memory allocation failed during parsing
    if (expr == NULL) return fail_msg("NULL expression.");

    switch (expr->type) {
    case E_SELECTION:
        assert(expr->value != NULL);
        selection->name = expr->value;
        return select_atoms(selection,expr->left,structure);
        break;
    case E_SYMBOL:
    case E_NAME:
    case E_RESN:
    case E_RESI:
    case E_CHAIN:
        return select_list(expr->type,selection,structure,expr->left);
        break;
    case E_AND:
    case E_OR: {
        struct selection *sl = selection_new(n), *sr = selection_new(n);
        if (sl && sr) {
            if ((ret = select_atoms(sl,expr->left,structure))) {
                if (ret == FREESASA_WARN) ++warn;
                if (ret == FREESASA_FAIL) ++err;
            }
            if ((ret = select_atoms(sr,expr->right,structure))) {
                if (ret == FREESASA_WARN) ++warn;
                if (ret == FREESASA_FAIL) ++err;
            }
            selection_join(selection,sl,sr,expr->type);
        } else {
            ++err;
        }
        selection_free(sl);
        selection_free(sr);
        if (err) return fail_msg("Error joining selections");
        break;
    }
    case E_NOT: {
        ret = select_atoms(selection,expr->right,structure);
        if (ret == FREESASA_WARN) ++warn;
        if (ret == FREESASA_FAIL) return FREESASA_FAIL;
        if (selection_not(selection)) return FREESASA_FAIL;
        break;
    }
    case E_ID:
    case E_NUMBER:
    case E_PLUS:
    case E_RANGE:
        // these four are handled by the RESN,SYMBOL,ETC
    default:
        return fail_msg("parser error");
    }
    if (warn) return FREESASA_WARN;
    return FREESASA_SUCCESS;
}

//for debugging
static void
print_expr(const expression *e,int level)
{
    if (freesasa_get_verbosity() == FREESASA_V_DEBUG) {
        fprintf(stderr,"\n");
        for (int i = 0; i < level; ++i) fprintf(stderr,"  ");
        if (e == NULL) fprintf(stderr,"()");
        else {
            fprintf(stderr,"(%s ",e_str(e->type));
            if (e->value) fprintf(stderr,": %s ",e->value);
            print_expr(e->left,level+1);
            print_expr(e->right,level+1);
            fprintf(stderr,")");
        }
        fflush(stderr);
    }
}

int
freesasa_select_area(const char *command,
                     char *name,
                     double *area,
                     const freesasa_structure *structure,
                     const freesasa_result *result)
{
    assert(name); assert(area); 
    assert(command); assert(structure); assert(result);
    assert(freesasa_structure_n(structure) == result->n_atoms);
    struct selection *selection = NULL;
    struct expression *expression = NULL;
    const int maxlen = FREESASA_MAX_SELECTION_NAME;
    double sasa = 0;
    int err = 0, warn = 0;
    *area = 0;
    name[0] = '\0';
    
    expression = get_expression(command);
    selection = selection_new(result->n_atoms);

    if (expression != NULL && selection != NULL) {
        switch (select_atoms(selection, expression, structure)) {
        case FREESASA_FAIL: 
            err = 1;
            break;
        case FREESASA_WARN: 
            warn = 1; // proceed with calculation, print warning later
        case FREESASA_SUCCESS: {
            for (int j = 0; j < selection->size; ++j) {
                sasa += selection->atom[j]*result->sasa[j];
            }
            
            *area = sasa;
            int len = strlen(selection->name);
            if (len > maxlen) {
                strncpy(name,selection->name,maxlen);
                name[maxlen+1] = '\0';
            }
            else {
                strcpy(name,selection->name);
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
        return freesasa_fail("in %s(): Problems parsing expression '%s'.",__func__,command);
    if (warn)
        return freesasa_warn("in %s(): There were warnings.",__func__);
    return FREESASA_SUCCESS;
}

int freesasa_selection_parse_error(expression *e,
                                   yyscan_t scanner,
                                   const char *msg)
{
    print_expr(e,0);
    if (freesasa_get_verbosity() == FREESASA_V_NORMAL) fprintf(stderr,"\n");
    return freesasa_fail(msg);
}
