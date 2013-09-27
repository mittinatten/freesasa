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

    //open PDB-file given as commandline argument
    FILE *pdb_file = fopen(argv[1],"r");

    //do the calculation using default parameters with PDB-file as input
    sasalib_calc_pdb(s,pdb_file);

    //print results
    sasalib_print(stdout,s);

    //clean up
    sasalib_free(s);
    fclose(pdb_file);

    return EXIT_SUCCESS;
}
