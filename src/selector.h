#ifndef SELECTOR_H
#define SELECTOR_H

typedef enum  
    {E_SELECTOR,E_SYMBOL,E_NAME,E_RESN,E_RESI,E_CHAIN,E_RANGE,E_ID,E_NUMBER}
expression_type;

typedef struct expression {
    struct expression *left;
    struct expression *right;
    expression_type type;
    const char *value;
} expression;

expression *
create_atom(expression_type type,
            const char* val);

expression *
create_selector(expression *e,
                const char* id);

expression *
create_selection(expression_type type,
                 expression *);

expression *
create_operation(expression_type type,
                 expression *left,
                 expression *right);

#endif /* SELECTOR_H */
