%{

#include "selection.h"
#include "parser.h"
#include "lexer.h"
    extern int freesasa_selection_parse_error(expression *e, yyscan_t scanner, const char *msg);
    int freesasa_yyerror(expression **expression, yyscan_t scanner, const char *msg) {
        return freesasa_selection_parse_error(*expression,scanner,msg);
    }

%}

%code requires {

#ifndef FREESASA_TYPEDEF_YY_SCANNER_T
#define FREESASA_TYPEDEF_YY_SCANNER_T
    typedef void* freesasa_yyscan_t;
#endif

}

%output "parser.c"
%defines "parser.h"
%name-prefix "freesasa_yy"
%define api.pure full
%lex-param { freesasa_yyscan_t scanner }
%parse-param {expression **expression }
%parse-param {freesasa_yyscan_t scanner }

%union {
    const char *value;
    expression *expression;
}

%token <value> T_NUMBER
%token <value> T_ID
%token <value> T_SELID

%token T_AND
%token T_OR
%token T_NOT

%token T_RESN
%token T_RESI
%token T_SYMBOL
%token T_NAME
%token T_CHAIN
%token T_MINUS

%precedence ATOM
%left T_OR
%left T_AND
%precedence T_NOT
%left '+'
%left '-'
%right T_MINUS

%type <expression> stmt
%type <expression> expr
%type <expression> list
%type <expression> r_range
%type <expression> c_range
%type <expression> id

%%

stmt:
  T_SELID ',' expr          { *expression = freesasa_selection_create($expr, $T_SELID); }
;

expr:
  '(' expr ')'           { $$ = $2; }
| expr T_AND expr        { $$ = freesasa_selection_operation(E_AND, $1, $3); }
| expr T_OR expr         { $$ = freesasa_selection_operation(E_OR, $1, $3); }
| T_NOT expr             { $$ = freesasa_selection_operation(E_NOT, NULL, $2); }
| T_RESN list            { $$ = freesasa_selection_selector(E_RESN, $list); }
| T_RESI r_range         { $$ = freesasa_selection_selector(E_RESI, $r_range); }
| T_SYMBOL list          { $$ = freesasa_selection_selector(E_SYMBOL, $list); }
| T_NAME list            { $$ = freesasa_selection_selector(E_NAME, $list); }
| T_CHAIN c_range        { $$ = freesasa_selection_selector(E_CHAIN, $c_range); }
;

list:
  id                     { $$ = $1; }
| id '+' list            { $$ = freesasa_selection_operation(E_PLUS, $1, $3); }
;

r_range:
  id                     { $$ = $1; }
| r_range '+' r_range    { $$ = freesasa_selection_operation(E_PLUS, $1, $3); }
| id '-' id              { $$ = freesasa_selection_operation(E_RANGE, $1, $3); }
| '-' id                 { $$ = freesasa_selection_operation(E_RANGE_OPEN_L, NULL, $2); }
| id '-'                 { $$ = freesasa_selection_operation(E_RANGE_OPEN_R, $1, NULL); }
;

c_range:
  id                     { $$ = $1; }
| c_range '+' c_range    { $$ = freesasa_selection_operation(E_PLUS, $1, $3); }
| id '-' id              { $$ = freesasa_selection_operation(E_RANGE, $1, $3); }
;

id:
  T_NUMBER               { $$ = freesasa_selection_atom(E_NUMBER, $1); }
| T_ID                   { $$ = freesasa_selection_atom(E_ID, $1); }
| T_MINUS T_NUMBER       { $$ = freesasa_selection_atom(E_NEGNUM, $2); }
;
