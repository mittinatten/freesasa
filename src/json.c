#if HAVE_CONFIG_H
  #include <config.h>
#endif
#include <json-c/json_object.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include "freesasa.h"
#include "freesasa_internal.h"
#include "classifier.h"

/** The functions in JSON-C does not seem to have any documented error
    return values. Therefore these errors are not caught. */

json_object *
freesasa_json_atom(freesasa_node *node,
                   int options)
{
    assert(node);
    json_object *atom = json_object_new_object();
    const freesasa_nodearea *area = freesasa_node_area(node);
    const char *name = freesasa_node_name(node);
    int n_len = strlen(name);
    char trim_name[n_len+1];

    sscanf(name, "%s", trim_name);

    json_object_object_add(atom, "name", json_object_new_string(trim_name));
    json_object_object_add(atom, "area", json_object_new_double(area->total));
    json_object_object_add(atom, "is-polar",
                           json_object_new_boolean(freesasa_node_atom_is_polar(node)));
    json_object_object_add(atom, "is-main-chain",
                           json_object_new_boolean(freesasa_atom_is_backbone(name)));
    json_object_object_add(atom, "radius",
                           json_object_new_double(freesasa_node_atom_radius(node)));

    return atom;
}

static json_object *
freesasa_json_nodearea(const freesasa_nodearea *area)
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
freesasa_json_residue(freesasa_node *node,
                      int options)
{
    assert(node);
    assert(freesasa_node_type(node) == FREESASA_NODE_RESIDUE);

    json_object *obj = json_object_new_object();
    const char *name = freesasa_node_name(node), *number;
    const freesasa_nodearea *abs = freesasa_node_area(node),
        *reference = freesasa_node_residue_reference(node);
    freesasa_nodearea rel;
    
    number = freesasa_node_residue_number(node);

    int n_len = strlen(number);
    char trim_number[n_len+1];
    sscanf(number, "%s", trim_number);

    json_object_object_add(obj, "name", json_object_new_string(name));
    json_object_object_add(obj, "number", json_object_new_string(trim_number));
    json_object_object_add(obj, "area", freesasa_json_nodearea(abs));

    if ((reference != NULL) && !(options & FREESASA_OUTPUT_SKIP_REL)) {
        freesasa_residue_rel_nodearea(&rel, abs, reference);
        json_object_object_add(obj, "relative-area", freesasa_json_nodearea(&rel));
    }
    
    json_object_object_add(obj, "n-atoms",
                           json_object_new_int(freesasa_node_residue_n_atoms(node)));
    return obj;
}

json_object *
freesasa_json_chain(freesasa_node *node,
                    int options)
{
    json_object *obj = json_object_new_object();
    const char *name = freesasa_node_name(node);

    json_object_object_add(obj, "label", json_object_new_string(name));
    json_object_object_add(obj, "n-residues", json_object_new_int(freesasa_node_chain_n_residues(node)));

    json_object_object_add(obj, "area",
                           freesasa_json_nodearea(freesasa_node_area(node)));

    return obj;
}

json_object *
freesasa_json_selection(const freesasa_selection **selections)
{
    assert(selections);
    json_object *obj = json_object_new_array();
    while (*selections) {
        json_object *json_selection = json_object_new_object();
        json_object_object_add(json_selection, "name", json_object_new_string(freesasa_selection_name(*selections)));
        json_object_object_add(json_selection, "area", json_object_new_double(freesasa_selection_area(*selections)));
        json_object_array_add(obj, json_selection);
        ++selections;
    };
    return obj;
}

json_object *
freesasa_json_structure(freesasa_node *node,
                        int options)
{
    json_object *obj = json_object_new_object();
    const freesasa_selection **selections = freesasa_node_structure_selections(node);
    
    json_object_object_add(obj, "chains", json_object_new_string(freesasa_node_structure_chain_labels(node)));
    json_object_object_add(obj, "model", json_object_new_int(freesasa_node_structure_model(node)));
    json_object_object_add(obj, "area",
                           freesasa_json_nodearea(freesasa_node_area(node)));
    if (selections != NULL) {
        json_object_object_add(obj, "selections",
                               freesasa_json_selection(selections));
    }
    return obj;
}

json_object *
freesasa_node2json(freesasa_node *node,
                   int exclude_type,
                   int options)
{
    json_object *obj, *array = NULL;
    int lowest = 0;
    int type = freesasa_node_type(node);
    freesasa_node *child = freesasa_node_children(node);
    if (child) {
        if (freesasa_node_type(child) == exclude_type) lowest = 1;
        if (!lowest) array = json_object_new_array();
    }

    switch (type) {
    case FREESASA_NODE_RESULT:
        obj = array;
        break;
    case FREESASA_NODE_STRUCTURE:
        obj = freesasa_json_structure(node, options);
        if (!lowest) json_object_object_add(obj, "chains", array);
        break;
    case FREESASA_NODE_CHAIN:
        obj = freesasa_json_chain(node, options);
        if (!lowest) json_object_object_add(obj, "residues", array);
        break;
    case FREESASA_NODE_RESIDUE:
        obj = freesasa_json_residue(node, options);
        if (!lowest) json_object_object_add(obj, "atoms", array);
        break;
    case FREESASA_NODE_ATOM:
        obj = freesasa_json_atom(node, options);
        break;
    case FREESASA_NODE_ROOT:
    default:
        assert(0 && "Tree illegal");
    }
    
    if (!lowest) {
        while (child) {
            json_object_array_add(array, freesasa_node2json(child, exclude_type, options));
            child = freesasa_node_next(child);
        }
    }
    return obj;
}

static json_object *
parameters2json(const freesasa_parameters *p)
{
    json_object *obj = json_object_new_object(), *res;

    json_object_object_add(obj, "algorithm", json_object_new_string(freesasa_alg_name(p->alg)));
    json_object_object_add(obj, "probe-radius", json_object_new_double(p->probe_radius));

    switch(p->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        res = json_object_new_int(p->shrake_rupley_n_points);
        break;
    case FREESASA_LEE_RICHARDS:
        res = json_object_new_int(p->lee_richards_n_slices);
        break;
    default:
        assert(0);
        break;
    }
    json_object_object_add(obj, "resolution", res);

    return obj;
}

static json_object *
json_result(freesasa_node *result,
            int options)
{
    json_object *obj = json_object_new_object();
    freesasa_nodetype exclude_type = FREESASA_NODE_NONE;
    const freesasa_parameters *parameters = freesasa_node_result_parameters(result);
    if (options & FREESASA_OUTPUT_STRUCTURE) exclude_type = FREESASA_NODE_CHAIN;
    if (options & FREESASA_OUTPUT_CHAIN) exclude_type = FREESASA_NODE_RESIDUE;
    if (options & FREESASA_OUTPUT_RESIDUE) exclude_type = FREESASA_NODE_ATOM;

    json_object_object_add(obj, "input", json_object_new_string(freesasa_node_name(result)));
    json_object_object_add(obj, "classifier", json_object_new_string(freesasa_node_classified_by(result)));
    json_object_object_add(obj, "parameters", parameters2json(parameters));
    json_object_object_add(obj, "structure", freesasa_node2json(result, exclude_type, options));

    return obj;
}

int
freesasa_write_json(FILE *output,
                    freesasa_node *root,
                    int options)

{
    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    json_object *results = json_object_new_array(),
        *json_root = json_object_new_object();
    freesasa_node *child = freesasa_node_children(root);

    json_object_object_add(json_root, "source", json_object_new_string(freesasa_string));
    json_object_object_add(json_root, "length-unit", json_object_new_string("Ångström"));
    json_object_object_add(json_root, "results", results);
    while (child) {
        json_object_array_add(results, json_result(child, options));
        child = freesasa_node_next(child);
    }

    fputs(json_object_to_json_string_ext(json_root, JSON_C_TO_STRING_PRETTY), output);
    json_object_put(json_root);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

