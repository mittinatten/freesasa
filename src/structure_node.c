#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "freesasa.h"
#include "freesasa_internal.h"
#include "classifier.h"

struct freesasa_structure_node {
    char *name;
    freesasa_node_type type;
    int first_atom, last_atom;
    const freesasa_structure *structure;
    freesasa_subarea *area;
    freesasa_structure_node *parent;
    freesasa_structure_node *children;
    freesasa_structure_node *next;
};

const freesasa_subarea freesasa_subarea_null = {NULL, 0, 0, 0, 0, 0};

static freesasa_structure_node *
structure_node_new(const char *name)
{
    freesasa_structure_node *node = malloc(sizeof(freesasa_structure_node));

    if (!node) goto memerr;

    *node = (freesasa_structure_node) {
        .name = NULL,
        .type = FREESASA_NODE_ATOM,
        .first_atom = 0, .last_atom = 0,
        .structure = NULL,
        .area = NULL,
        .parent = NULL,
        .children = NULL,
        .next = NULL
    };

    node->name = strdup(name);

    if (!node->name) goto memerr;
    
    return node;

 memerr:
    free(node);
    mem_fail();
    return NULL;
}

static void
structure_node_free(freesasa_structure_node *node)
{
    if (node) {
        freesasa_structure_node *current = node->children, *next;
        while (current) {
            next = current->next;
            structure_node_free(current);
            current = next;
        }
        free(node->name);
        free(node->area);
        free(node);
    }
}

typedef freesasa_structure_node* (*node_generator)(const freesasa_structure*,
                                                   int index);

static freesasa_structure_node *
structure_node_gen_children(freesasa_structure_node* parent,
                            int first,
                            int last,
                            node_generator ng)
{
    freesasa_structure_node *child, *first_child;

    first_child = ng(parent->structure, first);
    if (!first_child) return NULL;

    first_child->parent = parent;
    child = parent->children = first_child;

    for (int i = first+1; i <= last; ++i) {
        child->next = ng(parent->structure, i);
        if (!child->next) return NULL;
        child = child->next;
        child->parent = parent;
    }
    child->next = NULL;

    return first_child;
}

static freesasa_structure_node *
structure_atom_node(const freesasa_structure *structure,
                    int atom_index)
{
    freesasa_structure_node *atom = 
        structure_node_new(freesasa_structure_atom_name(structure, atom_index));
    
    if (!atom) {
        fail_msg("");
        return NULL;
    }
    atom->type = FREESASA_NODE_ATOM;
    atom->structure = structure;
    atom->first_atom = atom->last_atom = atom_index;

    return atom;
}

static freesasa_structure_node *
structure_residue_node(const freesasa_structure* structure,
                       int residue_index)
{
    freesasa_structure_node *residue = NULL;
    int first, last;

    residue = 
        structure_node_new(freesasa_structure_residue_name(structure, residue_index));
    
    if (!residue) {
        fail_msg("");
        return NULL;
    }

    freesasa_structure_residue_atoms(structure, residue_index, &first, &last);
    residue->type = FREESASA_NODE_RESIDUE;
    residue->structure = structure;
    residue->first_atom = first;
    residue->last_atom = last;
    
    if (!structure_node_gen_children(residue, first, last, structure_atom_node)) {
        structure_node_free(residue);
        return NULL;
    }
    
    return residue;
}

static freesasa_structure_node *
structure_chain_node(const freesasa_structure *structure,
                     int chain_index)
{
    const char *chains = freesasa_structure_chain_labels(structure);
    char name[2] = {chains[chain_index], '\0'};
    freesasa_structure_node *chain = NULL;
    int first_atom, last_atom, first_residue, last_residue;
    
    freesasa_structure_chain_atoms(structure, chains[chain_index], 
                                   &first_atom, &last_atom);

    chain = structure_node_new(name);
    if (!chain) {
        fail_msg("");
        return NULL;
    }

    chain->type = FREESASA_NODE_CHAIN;
    chain->structure = structure;
    chain->first_atom = first_atom;
    chain->last_atom = last_atom;

    freesasa_structure_chain_residues(structure, name[0],
                                      &first_residue, &last_residue);    

    if (!structure_node_gen_children(chain, first_residue, last_residue,
                                     structure_residue_node)) {
        structure_node_free(chain);
        return NULL;
    }

    return chain;
}

static int
structure_tree_add_areas(freesasa_structure_node *node,
                         const freesasa_result *result,
                         const freesasa_classifier *polar_classifier)
{
    assert(node);
    assert(result);

    if (polar_classifier == NULL) polar_classifier = &freesasa_default_classifier;
    freesasa_structure_node *child = node->children;
    freesasa_subarea atom,  *area = malloc(sizeof(freesasa_subarea));
    if (!area) return mem_fail();

    *area = freesasa_subarea_null;
    area->name = node->name;
    node->area = area;

    if (child) {
        while(child) {
            structure_tree_add_areas(child, result, polar_classifier);
            freesasa_add_subarea(node->area, child->area);
            child = child->next;
        }
    } else {
        for (int i = node->first_atom; i <= node->last_atom; ++i) {
            freesasa_atom_subarea(&atom, node->structure, result,
                                  polar_classifier, i);
            freesasa_add_subarea(node->area, &atom);
        }
    }
    return FREESASA_SUCCESS;
}

freesasa_structure_node *
freesasa_structure_tree(const freesasa_structure *structure,
                        const freesasa_result *result,
                        const freesasa_classifier *polar_classifier,
                        const char *name)
{
    freesasa_structure_node *root = structure_node_new(name);
    int n_chains = freesasa_structure_n_chains(structure);

    if (!root) goto cleanup;

    root->structure = structure;
    root->first_atom = 0;
    root->last_atom = freesasa_structure_n(structure);
    root->type = FREESASA_NODE_STRUCTURE;

    if (!structure_node_gen_children(root, 0, n_chains-1, structure_chain_node))
        goto cleanup;
    if (structure_tree_add_areas(root, result, polar_classifier))
        goto cleanup;

    return root;

 cleanup:
    structure_node_free(root);
    fail_msg("");
    return NULL;
}

int
freesasa_structure_tree_free(freesasa_structure_node *root) 
{
    if (root->parent) return fail_msg("Can't free node that isn't the root of its tree");
    
    structure_node_free(root);

    return FREESASA_SUCCESS;
}

const freesasa_subarea *
freesasa_structure_node_area(const freesasa_structure_node *node) 
{
    return node->area;
}

freesasa_structure_node *
freesasa_structure_node_children(freesasa_structure_node *node)
{
    return node->children;
}

freesasa_structure_node *
freesasa_structure_node_next(freesasa_structure_node *node)
{
    return node->next;
}

freesasa_structure_node *
freesasa_structure_node_parent(freesasa_structure_node *node)
{
    return node->parent;
}

freesasa_node_type
freesasa_structure_node_type(const freesasa_structure_node *node)
{
    return node->type;
}

const char *
freesasa_structure_node_name(const freesasa_structure_node *node)
{
    return node->name;
}

const freesasa_structure*
freesasa_structure_node_structure(const freesasa_structure_node *node)
{
    return node->structure;
}

void
freesasa_structure_node_atoms(const freesasa_structure_node *node,
                              int *first,
                              int *last)
{
    *first = node->first_atom;
    *last = node->last_atom;
}

void
freesasa_atom_subarea(freesasa_subarea *area,
                      const freesasa_structure *structure,
                      const freesasa_result *result,
                      const freesasa_classifier *polar_classifier,
                      int atom_index)
{
    const char *resn = freesasa_structure_atom_res_name(structure, atom_index);
    double a = result->sasa[atom_index];

    area->main_chain = area->side_chain = area->polar = area->apolar = 0;

    area->name = freesasa_structure_atom_name(structure, atom_index);
    area->total = a; 

    if (freesasa_atom_is_backbone(area->name))
        area->main_chain = a;
    else area->side_chain = a;

    if (freesasa_classifier_class(polar_classifier, resn, area->name))
        area->polar = a;
    else area->apolar = a;
}

void
freesasa_add_subarea(freesasa_subarea *sum,
                     const freesasa_subarea *term)
{
    sum->total += term->total;
    sum->side_chain += term->side_chain;
    sum->main_chain += term->main_chain;
    sum->polar += term->polar;
    sum->apolar += term->apolar;
}
