#ifndef FREESASA_JSON_H
#define FREESASA_JSON_H 1
#include <json-c/json_object.h>
#include "freesasa.h"

json_object *
freesasa_json_structure_tree(freesasa_structure_node *node,
                             const freesasa_classifier *classifier);

#endif
