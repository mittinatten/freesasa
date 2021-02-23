#ifndef CIF_HH
#define CIF_HH
#if HAVE_CONFIG_H
#  include <config.h>
#endif
#include <cstdio>
#include <vector>
#include "freesasa.h"

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options);

std::vector<freesasa_structure*>
freesasa_cif_structure_array(std::FILE *input,
                         int *n,
                         const freesasa_classifier *classifier,
                         int options);

// TODO this declaration does not have its implementation linked against. (Dont know why.)
int freesasa_warn(const char *format, ...);


#endif /* CIF_HH */