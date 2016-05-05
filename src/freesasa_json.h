#ifndef FREESASA_JSON_H
#define FREESASA_JSON_H 1
#include <json-c/json_object.h>
#include "freesasa.h"


json_object *
freesasa_json_atom(const freesasa_result *result,
                   const freesasa_structure *structure,
                   const freesasa_rsa_reference *rsa,
                   int atom_index);

json_object *
freesasa_json_residue(const freesasa_result *result,
                      const freesasa_structure *structure,
                      const freesasa_rsa_reference *rsa,
                      int residue_index,
                      freesasa_subarea *chain_area);
json_object *
freesasa_json_chain(const freesasa_result *result,
                    const freesasa_structure *structure,
                    const freesasa_rsa_reference *rsa,
                    char chain,
                    freesasa_subarea *structure_area);

json_object *
freesasa_json_result(const freesasa_result *result,
                     const freesasa_structure *structure,
                     const freesasa_rsa_reference *rsa,
                     const char *name);

#endif
