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

%precedence ATOM
%left T_OR
%left T_AND
%precedence T_NOT
%left '+'
%left '-'
%right UNARY

%type <expression> stmt
%type <expression> expr
%type <expression> list
%type <expression> range
%type <expression> atom

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
| T_RESI range           { $$ = freesasa_selection_selector(E_RESI, $range); }
| T_SYMBOL list          { $$ = freesasa_selection_selector(E_SYMBOL, $list); }
| T_NAME list            { $$ = freesasa_selection_selector(E_NAME, $list); }
| T_CHAIN range          { $$ = freesasa_selection_selector(E_CHAIN, $range); }
;

list:
  atom                   { $$ = $1; }
| atom '+' list          { $$ = freesasa_selection_operation(E_PLUS, $1, $3); }
;

range:
  atom                   { $$ = $1; }
| range '+' range        { $$ = freesasa_selection_operation(E_PLUS, $1, $3); }
| atom '-' atom          { $$ = freesasa_selection_operation(E_RANGE, $1, $3); }
;

atom:
  T_NUMBER               { $$ = freesasa_selection_atom(E_NUMBER, $1); }
| T_ID                   { $$ = freesasa_selection_atom(E_ID, $1); }
| '-' T_NUMBER %prec UNARY  { $$ = freesasa_selection_atom(E_NEGNUM, $2); }
;
