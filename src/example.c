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

#include <stdio.h>
#include <stdlib.h>

#include "freesasa.h"

/** \{ */
int main(int argc, char **argv)
{
    freesasa_structure *structure = NULL;
    freesasa_result *result = NULL;
    const freesasa_classifier *classifier = &freesasa_default_classifier;
    freesasa_nodearea area;

    /* Read structure from stdin */
    structure = freesasa_structure_from_pdb(stdin, classifier, 0);

    /* Calculate SASA using structure */
    if (structure) {
        result = freesasa_calc_structure(structure, NULL);
    }

    /* Calculate area of classes (Polar/Apolar/..) */
    if (result) {
        area = freesasa_result_classes(structure, result);
    } else {
        /* If there was an error at any step, we will end up here */
        printf("Error calculating SASA\n");
    }

    /* Print results */
    printf("Total  : %f A2\n", area.total);
    printf("Apolar : %f A2\n", area.apolar);
    printf("Polar  : %f A2\n", area.polar);

    /* Free allocated resources */
    freesasa_result_free(result);
    freesasa_structure_free(structure);

    return EXIT_SUCCESS;
}
/** \} */
