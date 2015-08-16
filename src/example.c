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

#include <stdlib.h>
#include "freesasa.h"

int main(int argc, char **argv) {

    freesasa_result result;
    freesasa_strvp *class_area;
    freesasa_structure *structure = freesasa_structure_from_pdb(stdin,0);
    double *radii = freesasa_structure_radius(structure,NULL);

    freesasa_calc_structure(&result,structure,radii,NULL);
    class_area = freesasa_result_classify(result,structure,NULL);

    printf("Total area: %f A2\n",result.total);
    for (int i = 0; i < class_area->n; ++i)
        printf("%s: %f A2\n",class_area->string[i],class_area->value[i]);

    // Free allocated resources, not strictly necessary but good style
    freesasa_result_free(result);
    freesasa_strvp_free(class_area);
    freesasa_structure_free(structure);
    free(radii);

    return EXIT_SUCCESS;
}
