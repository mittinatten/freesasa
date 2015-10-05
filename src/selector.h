#ifndef SELECTOR_H
#define SELECTOR_H

typedef enum  
    {E_SELECTOR,E_SYMBOL,E_NAME,E_RESN,E_RESI,E_CHAIN,E_ID,E_NUMBER,E_AND,E_OR,E_NOT,E_PLUS,E_RANGE}
expression_type;
static const char *e_str[] = 
    {"E_SELECTOR","E_SYMBOL","E_NAME","E_RESN","E_RESI","E_CHAIN",
     "E_ID","E_NUMBER","E_AND","E_OR","E_NOT","E_PLUS","E_RANGE"};

typedef struct expression {
    struct expression *left;
    struct expression *right;
    expression_type type;
    char *value;
} expression;

expression *
create_atom(expression_type type,
            const char* val);

expression *
create_selector(expression *selection,
                const char* id);

expression *
create_selection(expression_type type,
                 expression *list);

expression *
create_operation(expression_type type,
                 expression *left,
                 expression *right);

#endif /* SELECTOR_H */
