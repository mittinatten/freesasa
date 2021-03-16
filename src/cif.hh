#ifndef CIF_HH
#define CIF_HH

#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cstdio>
#include <vector>

#include "freesasa.h"
#include "freesasa_internal.h"

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options);

std::vector<freesasa_structure *>
freesasa_cif_structure_array(std::FILE *input,
                             int *n,
                             const freesasa_classifier *classifier,
                             int options);

int freesasa_write_cif(FILE *output,
                       freesasa_node *root,
                       int options);

#endif /* CIF_HH */