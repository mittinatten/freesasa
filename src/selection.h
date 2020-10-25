#ifndef SELECTION_H
#define SELECTION_H

typedef enum { E_SELECTION,
               E_SYMBOL,
               E_NAME,
               E_RESN,
               E_RESI,
               E_CHAIN,
               E_ID,
               E_NUMBER,
               E_NEGNUM,
               E_AND,
               E_OR,
               E_NOT,
               E_PLUS,
               E_RANGE,
               E_RANGE_OPEN_L,
               E_RANGE_OPEN_R } expression_type;

typedef struct expression {
    struct expression *left;
    struct expression *right;
    expression_type type;
    char *value;
} expression;

/** Create an atomic value E_ID or E_NUMBER */
expression *
freesasa_selection_atom(expression_type type,
                        const char *val);

/** Create a top level selection (E_SELECTION) */
expression *
freesasa_selection_create(expression *selection,
                          const char *id);

/** Create a property selector (E_SYMBOL, E_NAME, E_RESN, E_RESI or E_CHAIN) */
expression *
freesasa_selection_selector(expression_type type,
                            expression *list);

/** Create an operation (E_AND, E_OR, E_NOT, E_PLUS or E_RANGE) */
expression *
freesasa_selection_operation(expression_type type,
                             expression *left,
                             expression *right);

#endif /* SELECTION_H */
