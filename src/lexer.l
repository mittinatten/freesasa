%{

#include "selector.h"
#include "parser.h"
#include "util.h"

#include <stdio.h>

%}

%option outfile="lexer.c" header-file="lexer.h"

%option reentrant noyywrap never-interactive nounistd
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

NUMBER  \-*[[:digit:]]+
ID      [[:alpha:]][[:alnum:]]*
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
{NUMBER}  { yylval->value = strdup(yytext); return T_NUMBER; }
{ID}      { yylval->value = strdup(yytext); return T_ID; }




