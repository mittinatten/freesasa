/*
  Copyright Simon Mitternacht 2013-2016.

  This file is part of FreeSASA.

  FreeSASA is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FreeSASA is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FreeSASA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SELECTION_H
#define SELECTION_H

typedef enum  
    {E_SELECTION,E_SYMBOL,E_NAME,E_RESN,E_RESI,E_CHAIN,E_ID,E_NUMBER,E_AND,E_OR,E_NOT,E_PLUS,E_RANGE}
expression_type;

typedef struct expression {
    struct expression *left;
    struct expression *right;
    expression_type type;
    char *value;
} expression;

/** Create an atomic value E_ID or E_NUMBER */
expression *
freesasa_selection_atom(expression_type type,
                        const char* val);

/** Create a top level selection (E_SELECTION) */
expression *
freesasa_selection_create(expression *selection,
                          const char* id);

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
