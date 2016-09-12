#if HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "freesasa_internal.h"
#include "classifier.h"

struct atom_properties {
    int is_polar;
    int is_bb;
    double radius;
    char *pdb_line;
};

struct residue_properties {
    int n_atoms;
    char *number;
    freesasa_nodearea *reference;
};

struct chain_properties {
    int n_residues;
};

struct structure_properties {
    int n_chains;
    int n_atoms;
    int model;
    char *chain_labels;
    freesasa_result *result;
};

struct result_properties {
    char *classified_by;
    freesasa_parameters parameters;
    int n_structures;
};

struct freesasa_result_node {
    char *name;
    freesasa_node_type type;
    union properties {
        struct atom_properties atom;
        struct residue_properties residue;
        struct chain_properties chain;
        struct structure_properties structure;
        struct result_properties result;
    } properties;
    freesasa_nodearea *area;
    freesasa_result_node *parent;
    freesasa_result_node *children;
    freesasa_result_node *next;
};

const freesasa_nodearea freesasa_nodearea_null = {NULL, 0, 0, 0, 0, 0, 0};

static freesasa_result_node *
result_node_new(const char *name)
{
    freesasa_result_node *node = malloc(sizeof(freesasa_result_node));

    if (node == NULL) {
        goto memerr;
    }

    *node = (freesasa_result_node) {
        .name = NULL,
        .type = FREESASA_NODE_ATOM,
        .area = NULL,
        .parent = NULL,
        .children = NULL,
        .next = NULL
    };

    if (name) {
        node->name = strdup(name);
        if (node->name == NULL) {
            goto memerr;
        }
    }
    return node;

 memerr:
    free(node);
    mem_fail();
    return NULL;
}

static void
result_node_free(freesasa_result_node *node)
{
    if (node != NULL) {
        freesasa_result_node *current = node->children, *next;
        while (current) {
            next = current->next;
            result_node_free(current);
            current = next;
        }
        free(node->name);
        free(node->area);

        switch(node->type) {
        case FREESASA_NODE_ATOM:
            free(node->properties.atom.pdb_line);
            break;
        case FREESASA_NODE_RESIDUE:
            free(node->properties.residue.reference);
            free(node->properties.residue.number);
            break;
        case FREESASA_NODE_STRUCTURE:
            free(node->properties.structure.chain_labels);
            freesasa_result_free(node->properties.structure.result);
            break;
        case FREESASA_NODE_RESULT:
            free(node->properties.result.classified_by);
            break;
        default:
            break;
        }
        free(node);
    }
}

typedef freesasa_result_node* (*node_generator)(const freesasa_structure*,
                                                const freesasa_result*,
                                                int index);

static int
result_node_add_area(freesasa_result_node *node,
                     const freesasa_structure *structure,
                     const freesasa_result *result)
{
    freesasa_result_node *child = NULL;

    if (node->type == FREESASA_NODE_RESULT || node->type == FREESASA_NODE_ATOM) {
        return FREESASA_SUCCESS;
    }

    node->area = malloc(sizeof(freesasa_nodearea));
    if (node->area == NULL) {
        return mem_fail();
    }

    *node->area = freesasa_nodearea_null;
    node->area->name = node->name;

    child = node->children;
    while (child) {
        freesasa_add_nodearea(node->area, child->area);
        child = child->next;
    }

    return FREESASA_SUCCESS;
}

static freesasa_result_node *
result_node_gen_children(freesasa_result_node* parent,
                         const freesasa_structure *structure,
                         const freesasa_result *result,
                         int first,
                         int last,
                         node_generator ng)
{
    freesasa_result_node *child, *first_child;

    first_child = ng(structure, result, first);

    if (first_child == NULL){
        fail_msg("");
        return NULL;
    }
    
    first_child->parent = parent;
    child = parent->children = first_child;

    for (int i = first+1; i <= last; ++i) {
        child->next = ng(structure, result, i);
        if (child->next == NULL){
            fail_msg("");
            return NULL;
        }
        child = child->next;
        child->parent = parent;
    }
    child->next = NULL;

    result_node_add_area(parent, structure, result);
    
    return first_child;
}

static freesasa_result_node *
result_node_atom(const freesasa_structure *structure,
                 const freesasa_result *result,
                 int atom_index)
{
    freesasa_result_node *atom = 
        result_node_new(freesasa_structure_atom_name(structure, atom_index));
    const char *line;
    
    if (atom == NULL) {
        fail_msg("");
        return NULL;
    }

    atom->type = FREESASA_NODE_ATOM;
    atom->properties.atom = (struct atom_properties) {
        .is_polar = freesasa_structure_atom_class(structure, atom_index) == FREESASA_ATOM_POLAR,
        .is_bb = freesasa_atom_is_backbone(atom->name),
        .radius = freesasa_structure_atom_radius(structure, atom_index),
        .pdb_line = NULL,
    };

    line = freesasa_structure_atom_pdb_line(structure, atom_index);
    if (line != NULL) {
        atom->properties.atom.pdb_line = strdup(line);
        if (atom->properties.atom.pdb_line == NULL) {
            mem_fail();
            goto cleanup;
        }
    }
    
    atom->area = malloc(sizeof(freesasa_nodearea));
    if (atom->area == NULL) {
        mem_fail();
        goto cleanup;
    }

    atom->area->name = atom->name;
    freesasa_atom_nodearea(atom->area, structure, result, atom_index);

    return atom;

 cleanup:
    result_node_free(atom);
    return NULL;
}

static freesasa_result_node *
result_node_residue(const freesasa_structure *structure,
                    const freesasa_result *result,
                    int residue_index)
{
    freesasa_result_node *residue = NULL;
    const freesasa_nodearea *ref;
    int first, last;

    residue = result_node_new(freesasa_structure_residue_name(structure, residue_index));
    
    if (residue == NULL) {
        fail_msg("");
        return NULL;
    }

    residue->type = FREESASA_NODE_RESIDUE;

    freesasa_structure_residue_atoms(structure, residue_index, &first, &last);
    residue->properties.residue.n_atoms = last - first + 1;

    residue->properties.residue.number = strdup(freesasa_structure_residue_number(structure, residue_index));
    if (residue->properties.residue.number == NULL) {
        mem_fail();
        goto cleanup;
    }

    ref = freesasa_structure_residue_reference(structure, residue_index);
    if (ref != NULL) {
        residue->properties.residue.reference = malloc(sizeof(freesasa_nodearea));
        if (residue->properties.residue.reference == NULL) {
            mem_fail();
            goto cleanup;
        }
        // TODO copy name string too
        *residue->properties.residue.reference = *ref;
    }

    if (result_node_gen_children(residue, structure, result, first,
                                 last, result_node_atom) == NULL) {
        goto cleanup;
    }

    return residue;

 cleanup:
    result_node_free(residue);
    return NULL;
}

static freesasa_result_node *
result_node_chain(const freesasa_structure *structure,
                  const freesasa_result *result,
                  int chain_index)
{
    const char *chains = freesasa_structure_chain_labels(structure);
    char name[2] = {chains[chain_index], '\0'};
    freesasa_result_node *chain = NULL;
    int first_atom, last_atom, first_residue, last_residue;

    assert(strlen(chains) > chain_index);
    
    freesasa_structure_chain_atoms(structure, chains[chain_index], 
                                   &first_atom, &last_atom);

    chain = result_node_new(name);
    if (chain == NULL) {
        fail_msg("");
        return NULL;
    }

    chain->type = FREESASA_NODE_CHAIN;
    freesasa_structure_chain_residues(structure, name[0],
                                      &first_residue, &last_residue);    
    chain->properties.chain.n_residues = last_residue - first_residue + 1;
    
    if (result_node_gen_children(chain, structure, result,
                                    first_residue, last_residue,
                                    result_node_residue) == NULL) {
        fail_msg("");
        result_node_free(chain);
        return NULL;
    }

    return chain;
}

static freesasa_result_node *
result_node_structure(const freesasa_structure *structure,
                      const freesasa_result *result,
                      int dummy_index)
{
    freesasa_result_node *result_node = NULL;
    result_node = result_node_new(freesasa_structure_chain_labels(structure));
    
    if (result_node == NULL) {
        fail_msg("");
        return NULL;
    }

    result_node->type = FREESASA_NODE_STRUCTURE;
    result_node->properties.structure.n_chains = freesasa_structure_n_chains(structure);
    result_node->properties.structure.n_atoms = freesasa_structure_n(structure);
    result_node->properties.structure.result = NULL;
    result_node->properties.structure.chain_labels = strdup(freesasa_structure_chain_labels(structure));

    if (result_node->properties.structure.chain_labels == NULL) {
        fail_msg("");
        goto cleanup;
    }

    result_node->properties.structure.result = freesasa_result_clone(result);

    if (result_node->properties.structure.result == NULL) {
        fail_msg("");
        goto cleanup;
    }
    
    if (result_node_gen_children(result_node, structure, result, 0,
                                    freesasa_structure_n_chains(structure)-1,
                                    result_node_chain) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return result_node;
 cleanup:
    result_node_free(result_node);
    return NULL;
}

freesasa_result_node *
freesasa_result_tree_new(void)
{
    freesasa_result_node *tree = result_node_new(NULL);
    if (tree != NULL) {
        tree->type = FREESASA_NODE_ROOT;
    }
    return tree;
}

freesasa_result_node *
freesasa_result_tree_init(const freesasa_result *result,
                          const freesasa_structure *structure,
                          const char *name)
{
    freesasa_result_node *tree = result_node_new(NULL);

    tree->type = FREESASA_NODE_ROOT;

    if (tree == NULL) {
        fail_msg("");
    } else if (freesasa_result_tree_add_result(tree, result, structure, name) == FREESASA_FAIL) {
        fail_msg("");
        freesasa_result_node_free(tree);
        tree = NULL;
    }

    return tree;
}

int
freesasa_result_tree_add_result(freesasa_result_node *tree,
                                const freesasa_result *result,
                                const freesasa_structure *structure,
                                const char *name)
{
    freesasa_result_node *result_node = result_node_new(name);
    
    if (result_node == NULL) {
        goto cleanup;
    }

    result_node->type = FREESASA_NODE_RESULT;
    result_node->properties.result.n_structures = 1;
    result_node->properties.result.classified_by = strdup(freesasa_structure_classifier_name(structure));
    if (result_node->properties.result.classified_by == NULL) {
        mem_fail();
        goto cleanup;
    }
        
    if (result_node_gen_children(result_node, structure, result, 0, 0,
                                 result_node_structure) == NULL) {
        goto cleanup;
    }

    result_node->next = tree->children;
    tree->children = result_node;

    return FREESASA_SUCCESS;

 cleanup:
    result_node_free(result_node);
    fail_msg("");
    return FREESASA_FAIL;
}

int
freesasa_result_tree_join(freesasa_result_node *tree1,
                          freesasa_result_node **tree2)
{
    assert(tree1); assert(tree2); assert(*tree2);
    assert(tree1->type == FREESASA_NODE_ROOT);
    assert((*tree2)->type == FREESASA_NODE_ROOT);

    freesasa_result_node *child = tree1->children;
    if (child != NULL) {
        while (child->next) child = child->next;
        child->next = (*tree2)->children;
    } else {
        tree1->children = (*tree2)->children;
    }
    // tree1 takes over ownership, tree2 is invalidated.
    free(*tree2);
    *tree2 = NULL;

    return FREESASA_SUCCESS;
}

int
freesasa_result_node_free(freesasa_result_node *root) 
{
    if (root) {
        if (root->parent)
            return fail_msg("Can't free node that isn't the root of its tree");
        result_node_free(root);
    }
    return FREESASA_SUCCESS;
}

const freesasa_nodearea *
freesasa_result_node_area(const freesasa_result_node *node)
{
    assert(node->type != FREESASA_NODE_ROOT);
    return node->area;
}

const freesasa_result_node *
freesasa_result_node_children(const freesasa_result_node *node)
{
    return node->children;
}

const freesasa_result_node *
freesasa_result_node_next(const freesasa_result_node *node)
{
    return node->next;
}

const freesasa_result_node *
freesasa_result_node_parent(const freesasa_result_node *node)
{
    return node->parent;
}

freesasa_node_type
freesasa_result_node_type(const freesasa_result_node *node)
{
    return node->type;
}

const char *
freesasa_result_node_name(const freesasa_result_node *node)
{
    return node->name;
}

const char*
freesasa_result_node_classified_by(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_RESULT);
    return node->properties.result.classified_by;
}

int
freesasa_result_node_atom_is_polar(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.is_polar;
}

int
freesasa_result_node_atom_is_mainchain(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.is_bb;
}

double
freesasa_result_node_atom_radius(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.radius;
}

const char *
freesasa_result_node_atom_pdb_line(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.pdb_line;
}

int
freesasa_result_node_residue_n_atoms(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->properties.residue.n_atoms;
}

const char *
freesasa_result_node_residue_number(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->properties.residue.number;
}

const freesasa_nodearea *
freesasa_result_node_residue_reference(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->properties.residue.reference;
}

int
freesasa_result_node_chain_n_residues(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_CHAIN);
    return node->properties.chain.n_residues;
}

int
freesasa_result_node_structure_n_chains(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.n_chains;
}

int
freesasa_result_node_structure_n_atoms(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.n_atoms;
}

int
freesasa_result_node_structure_model(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.model;
}

const char *
freesasa_result_node_structure_chain_labels(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.chain_labels;
}

const freesasa_result *
freesasa_result_node_structure_result(const freesasa_result_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.result;
}

int
freesasa_atom_nodearea(freesasa_nodearea *area,
                       const freesasa_structure *structure,
                       const freesasa_result *result,
                       int atom_index)
{
    double a = result->sasa[atom_index];
    *area = freesasa_nodearea_null;

    area->total = a;

    if (freesasa_atom_is_backbone(freesasa_structure_atom_name(structure, atom_index)))
        area->main_chain = a;
    else area->side_chain = a;

    switch(freesasa_structure_atom_class(structure, atom_index)) {
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

    return FREESASA_SUCCESS;
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
                        int first_atom,
                        int last_atom)
{
    assert(area);
    assert(structure); assert(result);
    assert(first_atom <= last_atom);

    freesasa_nodearea term = freesasa_nodearea_null;
    for (int i = first_atom; i <= last_atom; ++i) {
        freesasa_atom_nodearea(&term, structure, result, i);
        freesasa_add_nodearea(area, &term);
    }
}
