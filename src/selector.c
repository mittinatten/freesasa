#include <string.h>
#include <stdlib.h>
#include "selector.h"
#include "parser.h"
#include "lexer.h"
#include "util.h"
#include "freesasa.h"

struct selection {
    int *atom;
    int size;
};

static expression *
expression_new() 
{
    struct expression *e = malloc(sizeof(expression));

    if (e == NULL) {mem_fail(); return NULL; }
    
    e->type = EMPTY;
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
create_atom(int type, const char* val)
{
    expression *e = expression_new();
    if (e == NULL) return NULL;

    e->type = E_ID;
    e->value = val;
}

expression *
create_selector(expresssion *e, const char* id)
{
    e->type = E_SELECTOR;
    e->value = id;
}

expression *
create_operation(int type, 
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
              const char *symbol)
{
    char atom[5];
    int n = strlen(symbol);
    for (int i = 0; i < selection->size; ++i) {
        sscanf(freesasa_structure_atom_name(structure,i),"%s",atom);
        if (strncmp(atom,symbol,n) == 0) selection->atom[i] = 1;
    }
}
static int
select_symbols(struct selection *selection,
               const freesasa_structure *structure,
               const char *symbols) 
{
    return FREESASA_SUCCESS;
}

static struct selection *
selection_join(const struct selection *s1,
               const struct selection *s2,
               int type)
{
    if (s1 == NULL || s2 == NULL) {
        freesasa_fail("%s: trying to join NULL selections\n",__func__);
        return NULL;
    }
    assert(s1->size == s2->size);
    int n = s1->size;
    struct selection *result = selection_new(n);
    if (result == NULL) {freesasa_fail(__func__); return NULL;}
    switch (type) {
    case AND:
        for (int i = 0; i < n; ++i) {
            result->atom[i] = s1->atom[i] && s2->atom[i];
        }
        break;
    case OR:
        for (int i = 0; i < n; ++i) {
            result->atom[i] = s1->atom[i] || s2->atom[i];
        }
        break;
    default: 
        assert(0);
    }

    return result;
}

// we want to just exit with empty arguments, since select_atoms() will then
// pass on a NULL pointer, error dealed with higher in hierarchy
static int
selection_not(struct selection *s)
{
    if (s == NULL) return FREESASA_FAIL;
    for (int i = 0; i < s->size; ++i) {
        s->atom[i] = ! s->atom[i];
    }
    return FREESASA_SUCCESS;
}

double
select_area(const char *selector,
            const freesasa_structure *structure,
            const freesasa_result *result)
{
    assert(selector); assert(structure); assert(result);
    assert(freesasa_structure_n(structure) == result->n_atoms);
    double sasa = 0;
    struct selection *selection;
    struct expression *expression = expression_new();
    if (expression == NULL) { 
        freesasa_fail(__func__); 
        return -1.0;
    }
    yyparse

    selection = select_atoms(expression,structure);
    if (selection == NULL) { 
        freesasa_fail(__func__); 
        return -1.0;
    }
    
    for (int i = 0; i < selection->size; ++i) {
        if (selection->atom[i]) sasa += result->sasa[i];
    }
    
    return sasa;
}


