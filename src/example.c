/*
  Copyright Simon Mitternacht 2013-2014.

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

#include <stdlib.h>
#include "freesasa.h"

int main(int argc, char **argv) {
    //initialize freesasa-object with default parameters
    freesasa_t *s = freesasa_init();

    //do the calculation using default parameters with protein
    //structure from STDIN
    freesasa_calc_pdb(s,stdin);

    //print results
    freesasa_log(stdout,s);
    printf("Total area: %f A2\n", freesasa_area_total(s));

    //clean up
    freesasa_free(s);

    return EXIT_SUCCESS;
}
