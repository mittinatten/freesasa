#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "freesasa.h"
#include "freesasa_internal.h"

typedef struct freesasa_structure_node freesasa_structure_node;

struct freesasa_structure_node {
    char *name;
    freesasa_node_type type;
    int first_atom, last_atom;
    const freesasa_structure *structure;
    freesasa_subarea *area;
    freesasa_structure_node **children;
    const freesasa_structure_node *parent;
};

const freesasa_subarea freesasa_subarea_null = {NULL, 0, 0, 0, 0, 0};

freesasa_structure_node *
structure_node_new(const char *name)
{
    const static freesasa_structure_node null_node = {
        .name = NULL,
        .type = FREESASA_NODE_ATOM,
        .first_atom = 0, .last_atom = 0,
        .structure = NULL,
        .area = NULL,
        .children = NULL,
        .parent = NULL
    };

    freesasa_structure_node *node = malloc(sizeof(freesasa_structure_node));

    if (!node) {
        mem_fail();
        return NULL;
    }

    *node = null_node;
    node->name = strdup(name);

    if (!node->name) {
        mem_fail();
        return NULL;
    }
    
    return node;
}

void
structure_node_free(freesasa_structure_node *node)
{
    if (node) {
        if (node->children) {
            for (int i = 0; node->children[i] != NULL; ++i) {
                structure_node_free(node->children[i]);
            }
        }
        free(node->children);
        free(node->name);
        free(node);
    }
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
structure_attach_atom_nodes(freesasa_structure_node *parent)
{
    int n_atoms = parent->last_atom - parent->first_atom + 1;
    freesasa_structure_node **atoms = 
        malloc(sizeof(freesasa_structure_node*) * (n_atoms + 1));

    if (!atoms) {
        fail_msg("");
        return NULL;
    }
    for (int i = 0; i < n_atoms + 1; ++i) atoms[i] = NULL;

    parent->children = atoms;
    for (int i = 0; i < n_atoms; ++i) {
        atoms[i] = structure_atom_node(parent->structure, i + parent->first_atom);
        if (!atoms[i]) return NULL;
        atoms[i]->parent = parent;
    }

    return parent;
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
    
    if (!structure_attach_atom_nodes(residue)) {
        structure_node_free(residue);
        return NULL;
    }
    
    return residue;
}

static freesasa_structure_node *
structure_attach_residue_nodes(freesasa_structure_node *parent,
                               char chain_label)
{
    int n_res, last, first;
    freesasa_structure_node **residues = NULL;
    
    freesasa_structure_chain_residues(parent->structure, chain_label,
                                      &first, &last);
    n_res = last - first + 1;
    residues = malloc(sizeof(freesasa_structure_node*) * (n_res+1));

    if (!residues) {
        mem_fail();
        return NULL;
    }
    for (int i = 0; i < n_res + 1; ++i) residues[i] = NULL;

    parent->children = residues;

    for (int i = 0; i < n_res; ++i) {
        residues[i] = structure_residue_node(parent->structure, i+first);
        if (!residues[i]) return NULL;
        residues[i]->parent = parent;
    }
    
    return parent;
}

static freesasa_structure_node *
structure_chain_node(const freesasa_structure *structure,
                     int chain_index)
{
    const char *chains = freesasa_structure_chain_labels(structure);
    char name[2] = {chains[chain_index], '\0'};
    freesasa_structure_node *chain = NULL;
    int first, last;
    
    freesasa_structure_chain_atoms(structure, chains[chain_index], &first, &last);

    chain = structure_node_new(name);
    if (!chain) {
        fail_msg("");
        return NULL;
    }

    chain->type = FREESASA_NODE_CHAIN;
    chain->structure = structure;
    chain->first_atom = first;
    chain->last_atom = last;

    if (!structure_attach_residue_nodes(chain, chains[chain_index])) {
        structure_node_free(chain);
        return NULL;
    }

    return chain;
}

static freesasa_structure_node *
structure_attach_chain_nodes(freesasa_structure_node *parent)
{
    int n_chains = freesasa_structure_n_chains(parent->structure);
    freesasa_structure_node **chains = 
        malloc(sizeof(freesasa_structure_node*) * (n_chains+1));

    if (!chains) {
        mem_fail();
        return NULL;
    }

    for (int i = 0; i < n_chains + 1; ++i) chains[i] = NULL;

    parent->children = chains;

    for (int i = 0; i < n_chains; ++i) {
        chains[i] = structure_chain_node(parent->structure, i);
        if (!chains[i]) return NULL;
        chains[i]->parent = parent;
    }

    return parent;
}

freesasa_structure_node *
freesasa_structure_tree_generate(freesasa_structure *structure,
                                 const char *name)
{
    freesasa_structure_node *root = structure_node_new(name);

    if (!root) {
        fail_msg("");
        return NULL;
    }

    root->structure = structure;
    root->first_atom = 0;
    root->last_atom = freesasa_structure_n(structure);
    root->type = FREESASA_NODE_STRUCTURE;

    if (!structure_attach_chain_nodes(root)) {
        structure_node_free(root);
        return NULL;
    }

    return root;
}

int
freesasa_structure_tree_free(freesasa_structure_node *root) 
{
    if (root->parent) return fail_msg("Can't free node that isn't the root of its tree");
    
    structure_node_free(root);

    return FREESASA_SUCCESS;
}

int
freesasa_structure_tree_fill(freesasa_structure_node *node,
                             const freesasa_result *result,
                             const freesasa_classifier *polar_classifier)
{
    if (polar_classifier == NULL) polar_classifier = &freesasa_default_classifier;
    freesasa_subarea atom,  *area = malloc(sizeof(freesasa_subarea));
    if (!area) return mem_fail();
    *area = freesasa_subarea_null;
    node->area = area;
    puts(node->name); fflush(stdout);
    if (node->children) {
        for (int i = 0; node->children[i] != NULL; ++i) {
            freesasa_structure_tree_fill(node->children[i], result, polar_classifier);
            freesasa_add_subarea(node->area, node->children[i]->area);
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

const freesasa_subarea *
freesasa_structure_node_area(const freesasa_structure_node *node) 
{
    return node->area;
}

const freesasa_structure_node const **
freesasa_structure_node_children(const freesasa_structure_node *node)
{
    return (const freesasa_structure_node **) node->children;
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

void
freesasa_atom_subarea(freesasa_subarea *area,
                      const freesasa_structure *structure,
                      const freesasa_result *result,
                      const freesasa_classifier *polar_classifier,
                      int atom_index)
{
    const char *resn = freesasa_structure_atom_res_name(structure, atom_index);
    const freesasa_classifier *bb = &freesasa_backbone_classifier;
    double a = result->sasa[atom_index];

    area->main_chain = area->side_chain = area->polar = area->apolar = 0;

    area->name = freesasa_structure_atom_name(structure, atom_index);
    area->total = a; 

    if (bb->sasa_class(resn, area->name, bb))
        area->main_chain = a;
    else area->side_chain = a;

    if (polar_classifier->sasa_class(resn, area->name, polar_classifier))
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
