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

#include <stdlib.h>
#include <stdio.h>
#include "freesasa.h"

int main(int argc, char **argv) {
    freesasa_result *result;
    freesasa_strvp *class_area;
    freesasa_structure *structure;

    /* Read structure from stdin */
    structure = freesasa_structure_from_pdb(stdin,NULL,0);

    /* Calculate SASA using structure and radii, store in
       'result'. NULL means default parameters. */
    result = freesasa_calc_structure(structure,NULL);
    
    /* Calculate area of classes (Polar/Apolar/..) using default
       classifier */
    class_area = freesasa_result_classify(result,structure,NULL);

    /* Print results */
    printf("Total area : %f A2\n",result->total);
    for (int i = 0; i < class_area->n; ++i)
        printf("%s : %f A2\n",class_area->string[i],
               class_area->value[i]);

    /* Free allocated resources, not strictly necessary in this
       context */
    freesasa_strvp_free(class_area);
    freesasa_result_free(result);
    freesasa_structure_free(structure);

    return EXIT_SUCCESS;
}
