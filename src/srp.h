/*
  Copyright Simon Mitternacht 2013-2015.

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

/**
   @file
   @author Simon Mitternacht
   
   This header provides arrays of test-points for the S&R
   algorithm. The arrays themselves are stored as hard-coded arrays
   and are referenced as const pointers. They were generated as
   described in the manual and are close to evenly distributed on the
   unit sphere. Sets of 20, 50, 100, 200, 500, 1000, 2000, and 5000
   points are availabe.
 */

/**
    Prints the legal values for the number of points as a
    comma-separated list, ending with newline. 

    This function is used to find out what are the available arrays of
    test points. Mainly intended for use by help functions writing to
    stdout. 

    @param output File to write the list to.
 */
void
freesasa_srp_print_n_opt(FILE *output);

/**
   Test if a given number of test-points is available.
   
   @param n Number of test points.
   @return 1 if n is allowed, 0 if not. 

 */
int
freesasa_srp_n_is_valid(int n);

/**
    Returns an array of n test points. 

    @param n Number of points.
    @return Array of coordinates of size 3*n, coordinates are stored as
    x1,y1,z1,x2,y2,z2,... Returns NULL if the argument n is invalid.
*/
const double*
freesasa_srp_get_points(int n);

#endif
