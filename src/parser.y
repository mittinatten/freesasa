%{

#include "selector.h"
#include "parser.h"
#include "lexer.h"

    int freesasa_yyerror(expression **expression, yyscan_t scanner, const char *msg) {
        return selector_parse_error(*expression,scanner,msg);
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

%type <expression> stmt
%type <expression> expr
%type <expression> list
%type <expression> range
%type <expression> atom

%%

stmt:
  T_ID ',' expr          { *expression = create_selector($expr,$T_ID); } 
;

expr:
  '(' expr ')'           { $$ = $2; }
| expr T_AND expr        { $$ = create_operation(E_AND, $1, $3); }
| expr T_OR expr         { $$ = create_operation(E_OR, $1, $3); }
| T_NOT expr             { $$ = create_operation(E_NOT, NULL, $2); }
| T_RESN list            { $$ = create_selection(E_RESN, $list); }
| T_RESI range           { $$ = create_selection(E_RESI, $range); }
| T_SYMBOL list          { $$ = create_selection(E_SYMBOL, $list); }
| T_NAME list            { $$ = create_selection(E_NAME, $list); }
| T_CHAIN range          { $$ = create_selection(E_CHAIN, $range); }
;

list:
  atom                   { $$ = $1; }
| atom '+' list          { $$ = create_operation(E_PLUS, $1, $3); }

range:
  atom                   { $$ = $1; }
| atom '-' atom          { $$ = create_operation(E_RANGE, $1, $3); }
| atom '-' atom '+' range{ $$ = create_operation(E_PLUS, create_operation(E_RANGE, $1, $3),$5); }
| atom '+' range         { $$ = create_operation(E_PLUS, $1, $3); }
;

atom:
  T_NUMBER               { $$ = create_atom(E_NUMBER,$1); }
| T_ID                   { $$ = create_atom(E_ID,$1); }
;
