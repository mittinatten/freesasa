#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "selector.h"
#include "parser.h"
#include "lexer.h"
#include "util.h"
#include "freesasa.h"

struct selection {
    const char* name;
    int *atom;
    int size;
};

extern freesasa_strvp* freesasa_strvp_new(int n);
extern void freesasa_strvp_free(freesasa_strvp *svp);
extern int yyparse(expression **expression, yyscan_t scanner);

static expression *
expression_new() 
{
    struct expression *e = malloc(sizeof(expression));

    if (e == NULL) {mem_fail(); return NULL; }
    
    e->type = E_SELECTOR;
    e->left = NULL;
    e->right = NULL;
    e->value = NULL;
   
    return e;
}

static void
expression_free(expression *expression)
{
    if (expression) {
        free(expression->left);
        free(expression->right);
        free(expression->value);
        free(expression);
    }
}

expression *
create_atom(expression_type type, const char* val)
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
    
    return e;
}

expression *
create_selector(expression *selection, const char* id)
{
    assert(id);
    expression *e = expression_new();
    if (e == NULL) return NULL;

    e->type = E_SELECTOR;
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
create_selection(expression_type type,
                 expression *list)
{
    expression *e = expression_new();
    if (e == NULL) return NULL;
    
    e->type = type;
    e->left = list;
    
    return e;
}

expression *
create_operation(expression_type type, 
                 expression *left, 
                 expression *right)
{
    expression *e = expression_new();
    if (e == NULL) return NULL;
    
    e->type = type;
    e->left = left;
    e->right = right;

    return e;
}


static expression *
get_expression(const char *selector) 
{
   yyscan_t scanner;
   YY_BUFFER_STATE state;
   int err;
   expression *expression = expression_new();
   if (yylex_init(&scanner) || expression == NULL) {
       freesasa_fail(__func__);
       return NULL;
   }
   state = yy_scan_string(selector, scanner);
   err = yyparse(&expression, scanner);
   if (err) {
       if (err == 1) freesasa_fail(__func__);
       if (err == 2) mem_fail();
       return NULL;
   }
   yy_delete_buffer(state, scanner);

   yylex_destroy(scanner);

   return expression;
}

static struct selection *
selection_new(int n)
{
    struct selection *selection = malloc(sizeof(struct selection));
    if (selection == NULL) { mem_fail(); return NULL; }
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

static void
select_symbol(struct selection *selection,
              const freesasa_structure *structure,
              const expression *expr)
{
    assert(expr);
    if (expr->type == E_PLUS) { 
        select_symbol(selection,structure,expr->right);
        select_symbol(selection,structure,expr->left);
    } else {
        assert(expr->type == E_ID);
        for (int i = 0; i < freesasa_structure_n(structure); ++i) {
            char atom[5];
            sscanf(freesasa_structure_atom_name(structure,i),"%s",atom);
            if (strncmp(atom,expr->value,strlen(expr->value)) == 0)
                selection->atom[i] = 1;
        }
    }
}

static int
selection_join(struct selection *target,
               const struct selection *s1,
               const struct selection *s2,
               int type)
{
    int n;
    if (s1 == NULL || s2 == NULL || target == NULL)
        return freesasa_fail("%s: trying to join NULL selections",__func__);
    
    assert(s1->size == s2->size);
    assert(s1->size == target->size);

    n = target->size;

    switch (type) {
    case E_AND:
        for (int i = 0; i < n; ++i) {
            target->atom[i] = s1->atom[i] && s2->atom[i];
        }
        break;
    case E_OR:
        for (int i = 0; i < n; ++i) {
            target->atom[i] = s1->atom[i] || s2->atom[i];
        }
        break;
    default: 
        assert(0);
    }

    return FREESASA_SUCCESS;
}

static int
selection_not(struct selection *s)
{
    if (s == NULL) return FREESASA_FAIL;
    for (int i = 0; i < s->size; ++i) {
        s->atom[i] = ! s->atom[i];
    }
    return FREESASA_SUCCESS;
}

/* Called recursively, the selection is built as we cover the expression tree */
int
select_atoms(struct selection* selection,
             const expression *expr,
             const freesasa_structure *structure)
{
    switch (expr->type) {
    case E_SELECTOR:
        assert(expr->value != NULL);
        selection->name = expr->value;
        select_atoms(selection,expr->left,structure);
        break;
    case E_SYMBOL:
        select_symbol(selection,structure,expr->left);
        break;
    case E_NAME:
        break;
    case E_RESN:
        break;
    case E_RESI:
        break;
    case E_CHAIN:
        break;
    case E_ID:
        break;
    case E_NUMBER:
        break;
    case E_AND:
    case E_OR: {
        int n = selection->size;
        struct selection *sl = selection_new(n),*sr = selection_new(n);
        if (sl && sr) {
            selection_join(selection,sl,sr,expr->type);
        } else {
            return freesasa_fail(__func__);
        }
        selection_free(sl);
        selection_free(sr);
        break;
    }
    case E_NOT:
        break;
    case E_PLUS:
    case E_RANGE:
        // these two are handled by the RESN,SYMBOL,ETC
    default:
        assert(0);
        return FREESASA_FAIL;
    }
    return FREESASA_SUCCESS;
}

freesasa_strvp*
select_area(const char **selector,
            int n_selector,
            const freesasa_structure *structure,
            const freesasa_result *result)
{
    assert(selector); assert(structure); assert(result);
    assert(freesasa_structure_n(structure) == result->n_atoms);
    double sasa = 0;
    struct selection *selection = NULL;
    struct expression *expression = NULL;
    freesasa_strvp *strvp = freesasa_strvp_new(n_selector);
    int err = 0;
    if (strvp == NULL) { mem_fail(); return NULL;}
   
    for (int i = 0; i < n_selector; ++i) {
        sasa = 0;
        assert(selector[i]);
        expression = get_expression(selector[i]);
        selection = selection_new(result->n_atoms);
        
        if (expression == NULL || selection == NULL) { ++err; break; }
        
        if (select_atoms(selection, expression, structure))
            { ++err; break; }
        
        assert(selection->atom); assert(selection->name);
        for (int j = 0; j < selection->size; ++j) {
            if (selection->atom[j] == 1) sasa += result->sasa[j];
        }
        strvp->value[i] = sasa;
        strvp->string[i] = strdup(selection->name);

        if (strvp->string[i] == NULL) {++err; break; }
        selection_free(selection);
        expression_free(expression);
    }
    if (err) {
        selection_free(selection);
        expression_free(expression);
        freesasa_strvp_free(strvp);
        return NULL;
    }
    return strvp;
}
