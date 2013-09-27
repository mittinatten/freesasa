/*
  Copyright Simon Mitternacht 2013.

  This file is part of Sasalib.
  
  Sasalib is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  Sasalib is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Sasalib.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "src/sasalib.h"

int main(int argc, char **argv) { 
    //initialize sasalib-object with default parameters
    sasalib_t *s = sasalib_init();

    //do the calculation using default parameters with protein
    //structure from STDIN
    sasalib_calc_pdb(s,stdin);

    //print results
    sasalib_log(stdout,s);
    printf("Total area: %f Å2\n", sasalib_total(s));

    //clean up
    sasalib_free(s);

    return EXIT_SUCCESS;
}
