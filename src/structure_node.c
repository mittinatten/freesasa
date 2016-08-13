#if HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "freesasa_internal.h"
#include "classifier.h"

struct freesasa_structure_node {
    char *name;
    char *classified_by; //! name of the classifier
    freesasa_node_type type;
    int first_atom, last_atom;
    const freesasa_structure *structure;
    freesasa_nodearea *area;
    freesasa_nodearea *reference; //! can be NULL
    freesasa_structure_node *parent;
    freesasa_structure_node *children;
    freesasa_structure_node *next;
};

const freesasa_nodearea freesasa_nodearea_null = {NULL, 0, 0, 0, 0, 0, 0};

static freesasa_structure_node *
structure_node_new(const char *name, const char *classified_by)
{
    freesasa_structure_node *node = malloc(sizeof(freesasa_structure_node));

    if (!node) goto memerr;

    *node = (freesasa_structure_node) {
        .name = NULL,
        .classified_by = NULL,
        .type = FREESASA_NODE_ATOM,
        .first_atom = 0, .last_atom = 0,
        .structure = NULL,
        .area = NULL,
        .reference = NULL,
        .parent = NULL,
        .children = NULL,
        .next = NULL
    };

    node->name = strdup(name);
    node->classified_by = strdup(classified_by);

    if (!node->name || !node->classified_by) goto memerr;
    
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
        free(node->classified_by);
        free(node->area);
        free(node->reference);
        free(node);
    }
}

typedef freesasa_structure_node* (*node_generator)(const freesasa_structure*,
                                                   const freesasa_classifier*,
                                                   int index);

static freesasa_structure_node *
structure_node_gen_children(freesasa_structure_node* parent,
                            const freesasa_structure *structure,
                            const freesasa_classifier *classifier,
                            int first,
                            int last,
                            node_generator ng)
{
    freesasa_structure_node *child, *first_child;

    first_child = ng(structure, classifier, first);
    if (!first_child) return NULL;

    first_child->parent = parent;
    child = parent->children = first_child;

    for (int i = first+1; i <= last; ++i) {
        child->next = ng(structure, classifier, i);
        if (!child->next) return NULL;
        child = child->next;
        child->parent = parent;
    }
    child->next = NULL;

    return first_child;
}

static freesasa_structure_node *
structure_atom_node(const freesasa_structure *structure,
                    const freesasa_classifier *classifier,
                    int atom_index)
{
    freesasa_structure_node *atom = 
        structure_node_new(freesasa_structure_atom_name(structure, atom_index),
                           freesasa_classifier_name(classifier));
    
    if (!atom) fail_msg("");
    else {
        atom->type = FREESASA_NODE_ATOM;
        atom->structure = structure;
        atom->first_atom = atom->last_atom = atom_index;
    }
    
    return atom;
}

static freesasa_structure_node *
structure_residue_node(const freesasa_structure* structure,
                       const freesasa_classifier* classifier,
                       int residue_index)
{
    freesasa_structure_node *residue = NULL;
    const freesasa_nodearea *ref;
    int first, last;

    residue = 
        structure_node_new(freesasa_structure_residue_name(structure, residue_index),
                           freesasa_classifier_name(classifier));
    
    if (!residue) {
        fail_msg("");
        return NULL;
    }

    freesasa_structure_residue_atoms(structure, residue_index, &first, &last);
    residue->type = FREESASA_NODE_RESIDUE;
    residue->structure = structure;
    residue->first_atom = first;
    residue->last_atom = last;
    ref = freesasa_classifier_residue_reference(classifier, residue->name);
    if (ref != NULL) {
        residue->reference = malloc(sizeof(freesasa_nodearea));
        if (residue->reference == NULL) {
            mem_fail();
            goto cleanup;
        }
        *residue->reference = *ref;
    }
    
    if (!structure_node_gen_children(residue, structure, classifier,
                                     first, last, structure_atom_node)) {
        goto cleanup;
    }
    
    return residue;

 cleanup:
    structure_node_free(residue);
    return NULL;
}

static freesasa_structure_node *
structure_chain_node(const freesasa_structure *structure,
                     const freesasa_classifier *classifier,
                     int chain_index)
{
    const char *chains = freesasa_structure_chain_labels(structure);
    char name[2] = {chains[chain_index], '\0'};
    freesasa_structure_node *chain = NULL;
    int first_atom, last_atom, first_residue, last_residue;
    
    freesasa_structure_chain_atoms(structure, chains[chain_index], 
                                   &first_atom, &last_atom);

    chain = structure_node_new(name, freesasa_classifier_name(classifier));
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

    if (!structure_node_gen_children(chain, structure, classifier,
                                     first_residue, last_residue,
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
    freesasa_nodearea atom,  *area = malloc(sizeof(freesasa_nodearea));
    if (!area) return mem_fail();

    *area = freesasa_nodearea_null;
    area->name = node->name;
    node->area = area;

    if (child) {
        while(child) {
            structure_tree_add_areas(child, result, polar_classifier);
            freesasa_add_nodearea(node->area, child->area);
            // store reference areas
            if (child->type == FREESASA_NODE_RESIDUE) {
                const freesasa_nodearea *ref =
                    freesasa_classifier_residue_reference(polar_classifier, node->name);
                if (ref != NULL) {
                    node->reference = malloc(sizeof(freesasa_nodearea));
                    if (node->reference) *node->reference = *ref;
                    else {
                        free(area);
                        return mem_fail();
                    }
                }
                
            }
            child = child->next;
        }
    } else {
        for (int i = node->first_atom; i <= node->last_atom; ++i) {
            freesasa_atom_nodearea(&atom, node->structure, result,
                                  polar_classifier, i);
            freesasa_add_nodearea(node->area, &atom);
        }
    }
    return FREESASA_SUCCESS;
}

freesasa_structure_node *
freesasa_result2tree(const freesasa_result *result,
                     const freesasa_structure *structure,
                     const freesasa_classifier *polar_classifier,
                     const char *name)
{
    if (polar_classifier == NULL) polar_classifier = &freesasa_default_classifier;
    freesasa_structure_node *root = structure_node_new(name, polar_classifier->name);
    int n_chains = freesasa_structure_n_chains(structure);

    if (!root) goto cleanup;

    root->structure = structure;
    root->first_atom = 0;
    root->last_atom = freesasa_structure_n(structure);
    root->type = FREESASA_NODE_STRUCTURE;

    if (!structure_node_gen_children(root, structure, polar_classifier,
                                     0, n_chains-1, structure_chain_node))
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
freesasa_structure_node_free(freesasa_structure_node *root) 
{
    if (root->parent) return fail_msg("Can't free node that isn't the root of its tree");
    
    structure_node_free(root);

    return FREESASA_SUCCESS;
}

const freesasa_nodearea *
freesasa_structure_node_area(const freesasa_structure_node *node)
{
    return node->area;
}

const freesasa_nodearea *
freesasa_structure_node_residue_reference(const freesasa_structure_node *node)
{
    return node->reference;
}

const freesasa_structure_node *
freesasa_structure_node_children(const freesasa_structure_node *node)
{
    return node->children;
}

const freesasa_structure_node *
freesasa_structure_node_next(const freesasa_structure_node *node)
{
    return node->next;
}

const freesasa_structure_node *
freesasa_structure_node_parent(const freesasa_structure_node *node)
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

const char*
freesasa_structure_node_classified_by(const freesasa_structure_node *node)
{
    return node->classified_by;
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
freesasa_atom_nodearea(freesasa_nodearea *area,
                       const freesasa_structure *structure,
                       const freesasa_result *result,
                       const freesasa_classifier *classifier,
                       int atom_index)
{
    const char *resn = freesasa_structure_atom_res_name(structure, atom_index);
    double a = result->sasa[atom_index];
    freesasa_atom_class atom_class;

    area->main_chain = area->side_chain = area->polar = area->apolar = 0;

    area->name = freesasa_structure_atom_name(structure, atom_index);
    area->total = a; 

    if (freesasa_atom_is_backbone(area->name))
        area->main_chain = a;
    else area->side_chain = a;

    switch(freesasa_classifier_class(classifier, resn, area->name)) {
    case FREESASA_ATOM_APOLAR:
        area->apolar = a;
        break;
    case FREESASA_ATOM_POLAR:
        area->polar = a;
        break;
    case FREESASA_ATOM_UNKNOWN:
        area->unknown = a;
        break;
    }
    
}

void
freesasa_add_nodearea(freesasa_nodearea *sum,
                      const freesasa_nodearea *term)
{
    sum->total += term->total;
    sum->side_chain += term->side_chain;
    sum->main_chain += term->main_chain;
    sum->polar += term->polar;
    sum->apolar += term->apolar;
    sum->unknown += term->unknown;
}

void
freesasa_range_nodearea(freesasa_nodearea *area,
                        const freesasa_structure *structure,
                        const freesasa_result *result,
                        const freesasa_classifier *classifier,
                        int first_atom,
                        int last_atom)
{
    assert(area);
    assert(structure); assert(result); assert(classifier);
    assert(first_atom <= last_atom);

    freesasa_nodearea term = freesasa_nodearea_null;
    for (int i = first_atom; i <= last_atom; ++i) {
        freesasa_atom_nodearea(&term, structure, result, classifier, i);
        freesasa_add_nodearea(area, &term);
    }
}
