%{

#include "selector.h"
#include "parser.h"
#include "util.h"

#include <stdio.h>

%}

%option outfile="lexer.c" header-file="lexer.h"
%option bison-bridge

LPAREN  "("
RPAREN  ")"
PLUS    "+"
DASH    "-"
COMMA   ","

AND     (?i:\&|and)
OR      (?i:\||or)
NOT     (?i:\!|not)

RESN    (?i:resn)
RESI    (?i:resi)
SYMBOL  (?i:symbol)
NAME    (?i:name)

NUMBER  \-*[0-9]+
ID      [a-zA-Z][a-zA-Z0-9]*
WS      [ \t\n\r]*

%%

{LPAREN}  { return T_LPAREN; }
{RPAREN}  { return T_RPAREN; }
{PLUS}    { return T_PLUS; }
{DASH}    { return T_DASH; }
{COMMA}   { return T_COMMA; }

{RESN}    { return T_RESN; }
{RESI}    { return T_RESI; }
{SYMBOL}  { return T_SYMBOL; }
{NAME}    { return T_NAME; }

{AND}     { return T_AND; }
{OR}      { return T_OR; }
{NOT}     { return T_NOT; }

{WS}      {}
{NUMBER}  { sscanf(yytext, "%d", &yylval->value); return T_NUMBER; }
{ID}      { sscanf(yytext, "%s", &yylval->value); return T_ID; }

%%

int yyerror(const char *msg) {
    return freesasa_fail(msg);
}

