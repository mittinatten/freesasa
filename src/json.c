#include <json-c/json_object.h>
#include <assert.h>
#include <string.h>
#include "freesasa.h"
#include "freesasa_internal.h"
#include "freesasa_json.h"

json_object *
freesasa_json_atom(const freesasa_result *result,
                   const freesasa_structure *structure,
                   const freesasa_rsa_reference *rsa,
                   int atom_index)
{
    json_object *atom = json_object_new_object();
    double radius = freesasa_structure_atom_radius(structure, atom_index);
    const char *aname = freesasa_structure_atom_name(structure, atom_index),
        *resn = freesasa_structure_atom_res_name(structure, atom_index);
    int name_len = strlen(aname);
    char trim_name[name_len + 1];
    int is_polar = rsa->polar_classifier->sasa_class(resn, aname, rsa->polar_classifier);
    int is_bb = rsa->bb_classifier->sasa_class(resn, aname, rsa->bb_classifier);
    double sasa = result->sasa[atom_index];
  
    sscanf(aname, "%s", trim_name); // get rid of whitespace
    
    json_object_object_add(atom, "name", json_object_new_string(trim_name));
    json_object_object_add(atom, "area", json_object_new_double(sasa));
    json_object_object_add(atom, "is-polar", json_object_new_boolean(is_polar));
    json_object_object_add(atom, "is-main-chain", json_object_new_boolean(is_bb));
    json_object_object_add(atom, "radius", json_object_new_double(radius));

    return atom;
}

static json_object *
freesasa_json_subarea(freesasa_subarea *area) 
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
freesasa_json_residue(const freesasa_result *result,
                      const freesasa_structure *structure,
                      const freesasa_rsa_reference *rsa,
                      int residue_index,
                      freesasa_subarea *chain_area)
{
    json_object *obj = json_object_new_object(), 
        *atom_array = json_object_new_array();
    freesasa_subarea abs, rel;
    const char *name = freesasa_structure_residue_name(structure, residue_index),
        *number = freesasa_structure_residue_number(structure, residue_index);
    int n_len = strlen(number);
    char trim_number[n_len+1];
    int first, last;

    sscanf(number, "%s", trim_number);
    json_object_object_add(obj, "name", json_object_new_string(name));
    json_object_object_add(obj, "number", json_object_new_string(trim_number));

    freesasa_rsa_val(&abs, &rel, residue_index, structure, result, rsa);
    freesasa_add_subarea(chain_area, &abs);
    json_object_object_add(obj, "abs", freesasa_json_subarea(&abs));
    if (rel.name)
        json_object_object_add(obj, "rel", freesasa_json_subarea(&rel));
    
    freesasa_structure_residue_atoms(structure, residue_index, &first, &last);
    json_object_object_add(obj, "n_atoms", json_object_new_int(last - first + 1));
    
    for (int i = first; i <= last; ++i) {
        json_object_array_add(atom_array, freesasa_json_atom(result, structure, rsa, i));
    }

    json_object_object_add(obj, "atoms", atom_array);
    return obj;
}

json_object *
freesasa_json_chain(const freesasa_result *result,
                    const freesasa_structure *structure,
                    const freesasa_rsa_reference *rsa,
                    char chain,
                    freesasa_subarea *structure_area)
{
    json_object *obj = json_object_new_object(),
        *residue_array = json_object_new_array();
    char label[2] = {chain, '\0'};
    freesasa_subarea chain_area = {label, 0, 0, 0, 0, 0};
    int first, last;

    json_object_object_add(obj, "label", json_object_new_string(label));
    freesasa_structure_chain_residues(structure, chain, &first, &last);
    json_object_object_add(obj, "n_residues", json_object_new_int(last - first + 1));

    for (int i = first; i <= last; ++i) {
        json_object_array_add(residue_array,
                              freesasa_json_residue(result, structure, rsa, i, &chain_area));
    }

    json_object_object_add(obj, "abs", freesasa_json_subarea(&chain_area));
    json_object_object_add(obj, "residues", residue_array);

    freesasa_add_subarea(structure_area, &chain_area);
    return obj;
}

json_object *
freesasa_json_result(const freesasa_result *result,
                     const freesasa_structure *structure,
                     const freesasa_rsa_reference *rsa,
                     const char *name)
{
    json_object *obj = json_object_new_object(),
        *chain_array = json_object_new_array();
    const char *chains = freesasa_structure_chain_labels(structure);
    int n_chains = strlen(chains);
    freesasa_subarea structure_area = {name, 0, 0, 0, 0, 0};
    
    json_object_object_add(obj, "name", json_object_new_string(name));
    json_object_object_add(obj, "n_chains", json_object_new_int(n_chains));
    for (int i = 0; i < n_chains; ++i) {
        json_object_array_add(chain_array,
                              freesasa_json_chain(result, structure, rsa, chains[i], &structure_area));
                              
    }
    json_object_object_add(obj, "abs", freesasa_json_subarea(&structure_area));

    json_object_object_add(obj, "chains", chain_array);

    return obj;
}
