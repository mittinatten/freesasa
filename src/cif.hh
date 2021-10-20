#ifndef CIF_HH
#define CIF_HH

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <vector>

#include "freesasa.h"
#include "freesasa_internal.h"

/**
    Generate ::freesasa_structure from CIF document.
    A copy of the document is stored in memory to allow reusing it later in
    ::freesasa_export_tree_to_cif.
    To free that memory call ::freesasa_cif_clear_docs();

    @param input Input file
    @param classifier Classifier to use
    @param structure_options Options (see ::freesasa_structure_add_atom)
*/
freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options);

/**
    Generate a set of structures from one file. See ::freesasa_structure_array.

    @param input Input file
    @param n Output variable for size
    @param classifier Classifier to use
    @param options Options (see ::freesasa_structure_array)
 */
std::vector<freesasa_structure *>
freesasa_cif_structure_array(std::FILE *input,
                             int *n,
                             const freesasa_classifier *classifier,
                             int options);

/**
    Output calculation results as part of a CIF document.
    Takes the original CIF document, and adds the data where appropriate.

    @param filename Output filename. If NULL output will be written to stdout.
    @param root Result tree to generate output from
*/
int freesasa_export_tree_to_cif(const char *filename,
                                freesasa_node *root);

#endif /* CIF_HH */
