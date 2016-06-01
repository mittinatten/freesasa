#if HAVE_CONFIG_H
  #include <config.h>
#endif
#include <json-c/json_object.h>
#include <assert.h>
#include <string.h>
#include "freesasa.h"
#include "freesasa_internal.h"
#include "classifier.h"

json_object *
freesasa_json_atom(const freesasa_structure_node *node)
{
    assert(node);
    json_object *atom = json_object_new_object();
    const freesasa_subarea *area = freesasa_structure_node_area(node);
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    int first, last;
    double radius;
    const char *name = freesasa_structure_node_name(node);
    int n_len = strlen(name);
    char trim_name[n_len+1];
    const char *resn = freesasa_structure_node_name(freesasa_structure_node_parent(node));
    int is_polar;
    int is_bb = freesasa_atom_is_backbone(name);
    double sasa = area->total;

    sscanf(name, "%s", trim_name);

    freesasa_structure_node_atoms(node, &first, &last);
    assert(first >= 0 && first < freesasa_structure_n(structure));
    radius = freesasa_structure_atom_radius(structure, first);
    is_polar = freesasa_structure_atom_class(structure, first);
  
    json_object_object_add(atom, "name", json_object_new_string(trim_name));
    json_object_object_add(atom, "area", json_object_new_double(sasa));
    json_object_object_add(atom, "is-polar", json_object_new_boolean(is_polar));
    json_object_object_add(atom, "is-main-chain", json_object_new_boolean(is_bb));
    json_object_object_add(atom, "radius", json_object_new_double(radius));

    return atom;
}

static json_object *
freesasa_json_subarea(const freesasa_subarea *area)
{
    json_object *obj = json_object_new_object();

    json_object_object_add(obj, "total",
                           json_object_new_double(area->total));
    json_object_object_add(obj, "polar",
                           json_object_new_double(area->polar));
    json_object_object_add(obj, "apolar",
                           json_object_new_double(area->apolar));
    json_object_object_add(obj, "main-chain",
                           json_object_new_double(area->main_chain));
    json_object_object_add(obj, "side-chain",
                           json_object_new_double(area->side_chain));

    return obj;
}

json_object *
freesasa_json_residue(const freesasa_structure_node *node)
{
    assert(node);
    assert(freesasa_structure_node_type(node) == FREESASA_NODE_RESIDUE);

    json_object *obj = json_object_new_object();
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node), *number;
    const freesasa_subarea *abs = freesasa_structure_node_area(node),
        *reference = freesasa_structure_node_residue_reference(node);
    freesasa_subarea rel;
    int first, last;

    freesasa_structure_node_atoms(node, &first, &last);
    number = freesasa_structure_atom_res_number(structure, first);

    int n_len = strlen(number);
    char trim_number[n_len+1];
    sscanf(number, "%s", trim_number);

    json_object_object_add(obj, "name", json_object_new_string(name));
    json_object_object_add(obj, "number", json_object_new_string(trim_number));
    json_object_object_add(obj, "abs", freesasa_json_subarea(abs));

    if (reference != NULL) {
        freesasa_residue_rel_subarea(&rel, abs, reference);
        json_object_object_add(obj, "rel", freesasa_json_subarea(&rel));
    }
    
    json_object_object_add(obj, "n_atoms", json_object_new_int(last - first + 1));
    return obj;
}

json_object *
freesasa_json_chain(const freesasa_structure_node *node)
{
    json_object *obj = json_object_new_object();
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node);
    int first, last;

    json_object_object_add(obj, "label", json_object_new_string(name));
    freesasa_structure_chain_residues(structure, name[0], &first, &last);
    json_object_object_add(obj, "n_residues", json_object_new_int(last - first + 1));

    json_object_object_add(obj, "abs",
                           freesasa_json_subarea(freesasa_structure_node_area(node)));

    return obj;
}

json_object *
freesasa_json_structure(const freesasa_structure_node *node)
{
    json_object *obj = json_object_new_object();
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node);
    int n_chains = freesasa_structure_n_chains(structure);

    json_object_object_add(obj, "name", json_object_new_string(name));
    json_object_object_add(obj, "n_chains", json_object_new_int(n_chains));
    json_object_object_add(obj, "abs",
                           freesasa_json_subarea(freesasa_structure_node_area(node)));

    return obj;
}

json_object *
freesasa_node2json(const freesasa_structure_node *node)
{
    json_object *obj, *array;
    const freesasa_structure_node *child = freesasa_structure_node_children(node);

    if (child) array = json_object_new_array();

    switch (freesasa_structure_node_type(node)) {
    case FREESASA_NODE_STRUCTURE:
        obj = freesasa_json_structure(node);
        json_object_object_add(obj, "chains", array);
        break;
    case FREESASA_NODE_CHAIN:
        obj = freesasa_json_chain(node);
        json_object_object_add(obj, "residues", array);
        break;
    case FREESASA_NODE_RESIDUE:
        obj = freesasa_json_residue(node);
        json_object_object_add(obj, "atoms", array);
        break;
    case FREESASA_NODE_ATOM:
        obj = freesasa_json_atom(node);
        break;
    default:
        assert(0 && "Tree illegal");
    }

    while (child) {
        json_object_array_add(array, freesasa_node2json(child));
        child = freesasa_structure_node_next(child);
    }

    return obj;
}


