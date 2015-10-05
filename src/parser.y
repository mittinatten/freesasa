%{

#include "selector.h"
#include "parser.h"
#include "lexer.h"
    
    int yyerror(expression **expression, yyscan_t scanner, const char *msg) {
        return selector_parse_error(*expression,scanner,msg);
    }

%}

%code requires {

#ifndef YY_TYPEDEF_YY_SCANNER_T
#define YY_TYPEDEF_YY_SCANNER_T
    typedef void* yyscan_t;
#endif

}

%output "parser.c"
%defines "parser.h"

%define api.pure full
%lex-param { yyscan_t scanner }
%parse-param {expression **expression }
%parse-param {yyscan_t scanner }

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
| T_RESI list            { $$ = create_selection(E_RESI, $list); }
| T_SYMBOL list          { $$ = create_selection(E_SYMBOL, $list); }
| T_NAME list            { $$ = create_selection(E_NAME, $list); }
;

list:
  atom                   { $$ = $1; }
| atom '-' atom          { $$ = create_operation(E_RANGE, $1, $3); }
| atom '-' atom '+' list { $$ = create_operation(E_PLUS, create_operation(E_RANGE, $1, $3),$5); }
| atom '+' list          { $$ = create_operation(E_PLUS, $1, $3); }
;

atom:
  T_NUMBER               { $$ = create_atom(E_NUMBER,$1); }
| T_ID                   { $$ = create_atom(E_ID,$1); }
;
