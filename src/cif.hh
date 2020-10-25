#ifndef CIF_HH
#define CIF_HH
#if HAVE_CONFIG_H
#  include <config.h>
#endif
#include <cstdio>
#include "freesasa.h"

freesasa_structure *
freesasa_structure_from_cif(std::FILE *input,
                            const freesasa_classifier *classifier,
                            int structure_options);

#endif /* CIF_HH */