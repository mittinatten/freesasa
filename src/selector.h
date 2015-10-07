#ifndef SELECTOR_H
#define SELECTOR_H

typedef enum  
    {E_SELECTOR,E_SYMBOL,E_NAME,E_RESN,E_RESI,E_CHAIN,E_ID,E_NUMBER,E_AND,E_OR,E_NOT,E_PLUS,E_RANGE}
expression_type;

static const char *e_str[] = 
    {"selector","symbol","name","resn","resi","chain",
     "id","number","and","or","not","plus","range"};

typedef struct expression {
    struct expression *left;
    struct expression *right;
    expression_type type;
    char *value;
} expression;

expression *
freesasa_selector_atom(expression_type type,
                       const char* val);

expression *
freesasa_selector_create(expression *selection,
                         const char* id);

expression *
freesasa_selector_selection(expression_type type,
                            expression *list);

expression *
freesasa_selector_operation(expression_type type,
                            expression *left,
                            expression *right);

#endif /* SELECTOR_H */
