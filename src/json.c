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
freesasa_json_atom(const freesasa_structure_node *node)
{
    assert(node);
    json_object *atom = json_object_new_object();
    const freesasa_nodearea *area = freesasa_structure_node_area(node);
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node);
    int first, last;
    int n_len = strlen(name);
    char trim_name[n_len+1];

    sscanf(name, "%s", trim_name);

    freesasa_structure_node_atoms(node, &first, &last);
    assert(first >= 0 && first < freesasa_structure_n(structure));
  
    json_object_object_add(atom, "name", json_object_new_string(trim_name));
    json_object_object_add(atom, "area", json_object_new_double(area->total));
    json_object_object_add(atom, "is-polar",
                           json_object_new_boolean(freesasa_structure_atom_class(structure, first) == FREESASA_ATOM_POLAR));
    json_object_object_add(atom, "is-main-chain",
                           json_object_new_boolean(freesasa_atom_is_backbone(name)));
    json_object_object_add(atom, "radius",
                           json_object_new_double(freesasa_structure_atom_radius(structure, first)));

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
freesasa_json_residue(const freesasa_structure_node *node)
{
    assert(node);
    assert(freesasa_structure_node_type(node) == FREESASA_NODE_RESIDUE);

    json_object *obj = json_object_new_object();
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node), *number;
    const freesasa_nodearea *abs = freesasa_structure_node_area(node),
        *reference = freesasa_structure_node_residue_reference(node);
    freesasa_nodearea rel;
    int first, last;

    freesasa_structure_node_atoms(node, &first, &last);
    number = freesasa_structure_atom_res_number(structure, first);

    int n_len = strlen(number);
    char trim_number[n_len+1];
    sscanf(number, "%s", trim_number);

    json_object_object_add(obj, "name", json_object_new_string(name));
    json_object_object_add(obj, "number", json_object_new_string(trim_number));
    json_object_object_add(obj, "area", freesasa_json_nodearea(abs));

    if (reference != NULL) {
        freesasa_residue_rel_nodearea(&rel, abs, reference);
        json_object_object_add(obj, "relative-area", freesasa_json_nodearea(&rel));
    }
    
    json_object_object_add(obj, "n-atoms", json_object_new_int(last - first + 1));
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
    json_object_object_add(obj, "n-residues", json_object_new_int(last - first + 1));

    json_object_object_add(obj, "area",
                           freesasa_json_nodearea(freesasa_structure_node_area(node)));

    return obj;
}

json_object *
freesasa_json_structure(const freesasa_structure_node *node)
{
    json_object *obj = json_object_new_object();
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    int n_chains = freesasa_structure_n_chains(structure);

    json_object_object_add(obj, "n-chains", json_object_new_int(n_chains));
    json_object_object_add(obj, "area",
                           freesasa_json_nodearea(freesasa_structure_node_area(node)));

    return obj;
}

json_object *
freesasa_node2json(const freesasa_structure_node *node, int exclude_type)
{
    json_object *obj, *array;
    int lowest = 0;
    int type = freesasa_structure_node_type(node);
    const freesasa_structure_node *child = freesasa_structure_node_children(node);
    if (child) {
        if (freesasa_structure_node_type(child) == exclude_type) lowest = 1;
        if (!lowest) array = json_object_new_array();
    }

    switch (type) {
    case FREESASA_NODE_ROOT:
        obj = json_object_new_array();
        break;
    case FREESASA_NODE_STRUCTURE:
        obj = freesasa_json_structure(node);
        if (!lowest) json_object_object_add(obj, "chains", array);
        break;
    case FREESASA_NODE_CHAIN:
        obj = freesasa_json_chain(node);
        if (!lowest) json_object_object_add(obj, "residues", array);
        break;
    case FREESASA_NODE_RESIDUE:
        obj = freesasa_json_residue(node);
        if (!lowest) json_object_object_add(obj, "atoms", array);
        break;
    case FREESASA_NODE_ATOM:
        obj = freesasa_json_atom(node);
        break;
    default:
        assert(0 && "Tree illegal");
    }
    
    if (!lowest) {
        while (child) {
            if (type != FREESASA_NODE_ROOT) {
                json_object_array_add(array, freesasa_node2json(child, exclude_type));
            } else {
                json_object_array_add(obj, freesasa_node2json(child, exclude_type));
            }
            child = freesasa_structure_node_next(child);
        }
    }
    return obj;
}

static json_object *
parameters2json(const freesasa_parameters *p)
{
    json_object *obj = json_object_new_object(), *res;
    extern const char *freesasa_alg_names[];

#ifdef HAVE_CONFIG_H
    json_object_object_add(obj, "source", json_object_new_string(PACKAGE_STRING));
#endif
    json_object_object_add(obj, "algorithm", json_object_new_string(freesasa_alg_names[p->alg]));
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

int
freesasa_write_json(FILE *output,
                    const freesasa_structure_node *root,
                    const freesasa_parameters *parameters,
                    int options)

{
    assert(freesasa_structure_node_type(root) == FREESASA_NODE_ROOT);

    json_object *obj = json_object_new_object(), *json_root = json_object_new_object();
    freesasa_node_type exclude_type = FREESASA_NODE_NONE;
    if (parameters == NULL) parameters = &freesasa_default_parameters;
    
    if (options & FREESASA_OUTPUT_STRUCTURE) exclude_type = FREESASA_NODE_CHAIN;
    if (options & FREESASA_OUTPUT_CHAIN) exclude_type = FREESASA_NODE_RESIDUE;
    if (options & FREESASA_OUTPUT_RESIDUE) exclude_type = FREESASA_NODE_ATOM;

    json_object_object_add(json_root, "FreeSASA-result", obj);

    json_object_object_add(obj, "input", json_object_new_string(freesasa_structure_node_name(root)));
    json_object_object_add(obj, "classifier", json_object_new_string(freesasa_structure_node_classified_by(root)));
    json_object_object_add(obj, "length-unit", json_object_new_string("Ångström"));
    json_object_object_add(obj, "parameters", parameters2json(parameters));
    json_object_object_add(obj, "structures", freesasa_node2json(root, exclude_type));

    fputs(json_object_to_json_string_ext(json_root, JSON_C_TO_STRING_PRETTY), output);
    json_object_put(json_root);

    fflush(output);
    if (ferror(output)) {
        return fail_msg(strerror(errno));
    }
    return FREESASA_SUCCESS;
}

