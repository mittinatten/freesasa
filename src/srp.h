/*
  Copyright Simon Mitternacht 2013.

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

#ifndef FREESASA_SRP_H
#define FREESASA_SRP_H

/** Prints the legal values for the number of points as a
    comma-separated list, ending with newline. */
void freesasa_srp_print_n_opt(FILE*);

/** Returns 1 if n-value is allowed, 0 if not. */
int freesasa_srp_n_is_valid(int n);

/** Returns an array of n test points (array has size 3*n). If n is
    not one of the legal values, an error message is printed and
    exit(1) is called. */
const double* freesasa_srp_get_points(int n);

#endif
