#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <assert.h>
#include "selector.h"
#include "parser.h"
#include "lexer.h"
#include "util.h"
#include "freesasa.h"
#include "pdb.h"
#include "classify.h"

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

    for (int i = 0; i < strlen(val); ++i) e->value[i] = toupper(val[i]);
    
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
    if (strncmp(resn, id, strlen(id)) == 0)
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
        default:
            assert(0);
            break;
        }
        if (match) selection->atom[i] = 1;
    }
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
    default:
        assert(0);
        break;
    }
    return FREESASA_SUCCESS;
}
static int
select_list(expression_type parent_type,
            struct selection *selection,
            const freesasa_structure *structure,
            const expression *expr)
{
    assert(expr);
    int lower, upper, resr, resl;
    expression *left = expr->left, *right = expr->right;
    switch(expr->type) {
    case E_PLUS: 
        resl = select_list(parent_type,selection,structure,left);
        resr = select_list(parent_type,selection,structure,right);
        if (resl == FREESASA_WARN || resr == FREESASA_WARN)
            return FREESASA_WARN;
        break;
    case E_RANGE:
        assert(parent_type == E_RESI);
        if (left->type != E_NUMBER || right->type != E_NUMBER) {
            return freesasa_warn("select: Range '%s-%s' invalid, needs to be two numbers. "
                                 "Will be ignored.",left->value,right->value);
        }
        lower = atoi(left->value);
        upper = atoi(right->value);
        for (int i = 0; i < selection->size; ++i) {
            int resi = atoi(freesasa_structure_atom_res_number(structure,i));
            if (resi >= lower && resi <= upper) 
                selection->atom[i] = 1;
        }
        break;
    case E_ID: case E_NUMBER:
        if (is_valid_id(parent_type,expr) == FREESASA_SUCCESS)
            select_id(parent_type,selection,structure,expr->value);
        break;
    default:
        freesasa_fail("%s: %s %s",__func__,e_str[parent_type],e_str[expr->type]);
        assert(0);
        break;
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
        return freesasa_fail("%s: trying to join NULL selections",__func__);
    
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
    if (s == NULL) return FREESASA_FAIL;
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
    int warn = 0;
    switch (expr->type) {
    case E_SELECTOR:
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
        int n = selection->size;
        struct selection *sl = selection_new(n),*sr = selection_new(n);
        if (sl && sr) {
            if (select_atoms(sl,expr->left,structure)  == FREESASA_WARN) ++warn;
            if (select_atoms(sl,expr->right,structure) == FREESASA_WARN) ++warn;
            selection_join(selection,sl,sr,expr->type);
        } else {
            return freesasa_fail(__func__);
        }
        selection_free(sl);
        selection_free(sr);
        break;
    }
    case E_NOT:
        if (select_atoms(selection,expr->left,structure) == FREESASA_WARN) ++warn;
        selection_not(selection);
        break;
    case E_ID:
    case E_NUMBER:
    case E_PLUS:
    case E_RANGE:
        // these four are handled by the RESN,SYMBOL,ETC
    default:
        return freesasa_fail("%s: parser error.",__func__);
    }
    if (warn) return FREESASA_WARN;
    return FREESASA_SUCCESS;
}

//for debugging
static void
print_expr(const expression *e,int level)
{
    fprintf(stderr,"\n");
    for (int i = 0; i < level; ++i) fprintf(stderr,"  ");
    if (e == NULL) fprintf(stderr,"()");
    else {
        fprintf(stderr,"(%s ",e_str[e->type]);
        if (e->value) fprintf(stderr,": %s ",e->value);
        print_expr(e->left,level+1);
        print_expr(e->right,level+1);
        fprintf(stderr,")");
    }
    fflush(stderr);
}

int
freesasa_select_area(const char *command,
                     char **name,
                     double *area,
                     const freesasa_structure *structure,
                     const freesasa_result *result)
{
    assert(name); assert(area); 
    assert(command); assert(structure); assert(result);
    assert(freesasa_structure_n(structure) == result->n_atoms);
    struct selection *selection = NULL;
    struct expression *expression = NULL;
    double sasa = 0;
    int err = 0, warn = 0, flag;
    *area = 0;
    *name = "";

    expression = get_expression(command);
    selection = selection_new(result->n_atoms);
    //print_expr(expression,0);
    if (expression != NULL && selection != NULL) {
        flag = select_atoms(selection, expression, structure);
        if (flag == FREESASA_FAIL) err = 1;
        if (flag == FREESASA_WARN) warn = 1;
        if (flag == FREESASA_SUCCESS || flag == FREESASA_WARN) {
            //int count = 0;
            for (int j = 0; j < selection->size; ++j) {
                sasa += selection->atom[j]*result->sasa[j];
                //count += selection->atom[j];
            }
            
            *area = sasa;
            *name = strdup(selection->name);
            //printf(">> %s %f %d\n",*name,*area,count); fflush(stdout);
            
            if (*name == NULL) {
                err = 1; 
                mem_fail();
            }
        }
    } else {
        err = 1;
    }
    selection_free(selection);
    expression_free(expression);
    
    if (err) return freesasa_fail(__func__);
    if (warn) return freesasa_warn("%s: there were warnings.",__func__);
    return FREESASA_SUCCESS;
}

int selector_parse_error(expression *e,
                         yyscan_t scanner,
                         const char *msg)
{
    print_expr(e,0);
    fprintf(stderr,"\n");
    return freesasa_fail("%s: %s: %s %s",__func__,msg,e_str[e->type],e->value);
}


