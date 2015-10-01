%{

#include "util.h"
#include "selector.h"
#include "parser.h"
#include "lexer.h"
    
    int yyerror(expression **expression, yyscan_t scanner, const char *msg) {
        return freesasa_fail(msg);
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

%token T_LPAREN
%token T_RPAREN
%token T_PLUS
%token T_DASH
%token T_COMMA

%token T_AND
%token T_OR
%token T_NOT

%token <value> T_NUMBER
%token <value> T_ID

%token T_RESN
%token T_RESI
%token T_SYMBOL
%token T_NAME

%left T_DASH
%left T_PLUS
%left T_AND
%left T_OR
%right T_NOT

%type <expression> selector
%type <expression> selection
%type <expression> list
%type <expression> i_list
%type <expression> atom

%%

input: selector ;

selector:
T_ID[ID] T_COMMA selection[S]  { *expression = create_selector($S,$ID); } 

list:
  atom                      { $$ = $1; } 
| atom[L] T_PLUS list[R]    { $$ = create_operation(E_PLUS, $L, $R); }
;

i_list:
  atom                      { $$ = $1; }  
| atom[L] T_DASH i_list[R]  { $$ = create_operation(E_RANGE, $L, $R); }
| atom[L] T_PLUS i_list[R]  { $$ = create_operation(E_PLUS, $L, $R); }
;

selection: 
  T_RESN list[L]            { $$ = create_selection(E_RESN, $L); }
| T_RESI i_list[L]          { $$ = create_selection(E_RESI, $L); }
| T_SYMBOL list[L]          { $$ = create_selection(E_SYMBOL, $L); }
| T_NAME list[L]            { $$ = create_selection(E_NAME, $L); }
| T_LPAREN selection[S] T_RPAREN 
                            { $$ = $S; }
| selection[L] T_AND selection[R]
                            { $$ = create_operation(E_AND, $L, $R); }
| selection[L] T_OR selection[R]
                            { $$ = create_operation(E_OR, $L, $R); }
| T_NOT selection[R]        { $$ = create_operation(E_NOT, NULL, $R); }
;

atom:
  T_ID                      { $$ = create_atom(E_ID,$1); }
| T_NUMBER                  { $$ = create_atom(E_NUMBER,$1); }
;
