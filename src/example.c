/**
    @file
    @author Simon Mitternacht
    @copyright [MIT License](md_license.html)

    @brief Short program that illustrates how to use the most basic
    functionality of the API

    The program does basic error handling, printing an unspecific
    error message if anything fails, in addition to the library
    errors.
 */

#include <stdlib.h>
#include <stdio.h>
#include "freesasa.h"

//! \{
int main(int argc, char **argv) {
    freesasa_structure *structure = NULL;
    freesasa_result *result = NULL;
    freesasa_strvp *class_area = NULL;

    /* Read structure from stdin */
    structure = freesasa_structure_from_pdb(stdin,NULL,0);

    /* Calculate SASA using structure */
    if (structure) {
        result = freesasa_calc_structure(structure,NULL);
    }

    /* Calculate area of classes (Polar/Apolar/..) */
    if (result) {
        class_area = freesasa_result_classify(result,structure,NULL);
    }
    
    /* Print results */
    if (class_area) {
        printf("Total area : %f A2\n",result->total);
        for (int i = 0; i < class_area->n; ++i)
            printf("%s : %f A2\n",class_area->string[i],
                   class_area->value[i]);
    } else {
        /* If there was an error at any step, we will end up here */
        printf("Error calculating SASA\n");
    }

    /* Free allocated resources */
    freesasa_strvp_free(class_area);
    freesasa_result_free(result);
    freesasa_structure_free(structure);

    return EXIT_SUCCESS;
}
//! \}
