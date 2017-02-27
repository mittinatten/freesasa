%{

#include "selection.h"
#include "parser.h"
#include "freesasa_internal.h"

#include <stdio.h>

%}

%option outfile="lexer.c" header-file="lexer.h"
%option prefix="freesasa_yy"
%option reentrant noyywrap never-interactive nounistd
%option bison-bridge

AND     (?i:\&|and)
OR      (?i:\||or)
NOT     (?i:\!|not)
MINUS   \\-

RESN    (?i:resn)
RESI    (?i:resi)
SYMBOL  (?i:symbol)
NAME    (?i:name)
CHAIN   (?i:chain)

NUMBER  [[:digit:]]+
ID      [[:alnum:]\_]+\'*
SELID   [[:alnum:]\_\-\+]+
WS      [ \t\n\r]*

%%

","       { return ','; }
"-"       { return '-'; }
"+"       { return '+'; }
"("       { return '('; }
")"       { return ')'; }

{RESN}    { return T_RESN; }
{RESI}    { return T_RESI; }
{SYMBOL}  { return T_SYMBOL; }
{NAME}    { return T_NAME; }
{CHAIN}   { return T_CHAIN; }

{AND}     { return T_AND; }
{OR}      { return T_OR; }
{NOT}     { return T_NOT; }
{MINUS}   { return T_MINUS; }

{WS}          {}
{NUMBER}      { yylval->value = strdup(yytext); return T_NUMBER; }
{ID}          { yylval->value = strdup(yytext); return T_ID; }
{SELID}/","   { yylval->value = strdup(yytext); return T_SELID; }
