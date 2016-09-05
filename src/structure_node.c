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
};

struct residue_properties {
    freesasa_nodearea *reference;
    char *number;
    int n_atoms;
};

struct chain_properties {
    int n_residues;
};

struct structure_properties {
    int n_chains;
    char *chain_labels;
};

struct tree_properties {
    char *classified_by;
    int n_structures;
};

struct freesasa_structure_node {
    char *name;
    freesasa_node_type type;
    union {
        struct atom_properties atom_prop;
        struct residue_properties residue_prop;
        struct chain_properties chain_prop;
        struct structure_properties structure_prop;
        struct tree_properties tree_prop;
    };
    freesasa_nodearea *area;
    freesasa_structure_node *parent;
    freesasa_structure_node *children;
    freesasa_structure_node *next;
};

const freesasa_nodearea freesasa_nodearea_null = {NULL, 0, 0, 0, 0, 0, 0};

static freesasa_structure_node *
structure_node_new(const char *name)
{
    freesasa_structure_node *node = malloc(sizeof(freesasa_structure_node));

    if (node == NULL) {
        goto memerr;
    }

    *node = (freesasa_structure_node) {
        .name = NULL,
        .type = FREESASA_NODE_ATOM,
        .area = NULL,
        .parent = NULL,
        .children = NULL,
        .next = NULL
    };

    node->name = strdup(name);
  
    if (node->name == NULL) {
        goto memerr;
    }
    return node;

 memerr:
    free(node);
    mem_fail();
    return NULL;
}

static void
structure_node_free(freesasa_structure_node *node)
{
    if (node != NULL) {
        freesasa_structure_node *current = node->children, *next;
        while (current) {
            next = current->next;
            structure_node_free(current);
            current = next;
        }
        free(node->name);
        free(node->area);
        switch(node->type) {
        case FREESASA_NODE_RESIDUE:
            free(node->residue_prop.reference);
            free(node->residue_prop.number);
            break;
        case FREESASA_NODE_ROOT:
            free(node->tree_prop.classified_by);
            break;
        case FREESASA_NODE_STRUCTURE:
            free(node->structure_prop.chain_labels);
            break;
        default:
            break;
        }
        free(node);
    }
}

typedef freesasa_structure_node* (*node_generator)(const freesasa_structure*,
                                                   const freesasa_classifier*,
                                                   const freesasa_result*,
                                                   int index);

static int
structure_node_add_area(freesasa_structure_node *node,
                        const freesasa_structure *structure,
                        const freesasa_classifier *classifier,
                        const freesasa_result *result)
{
    freesasa_structure_node *child = NULL;

    if (node->type == FREESASA_NODE_ROOT || node->type == FREESASA_NODE_ATOM) {
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

static freesasa_structure_node *
structure_node_gen_children(freesasa_structure_node* parent,
                            const freesasa_structure *structure,
                            const freesasa_classifier *classifier,
                            const freesasa_result *result,
                            int first,
                            int last,
                            node_generator ng)
{
    freesasa_structure_node *child, *first_child;

    first_child = ng(structure, classifier, result, first);

    if (first_child == NULL){
        fail_msg("");
        return NULL;
    }
    
    first_child->parent = parent;
    child = parent->children = first_child;

    for (int i = first+1; i <= last; ++i) {
        child->next = ng(structure, classifier, result, i);
        if (child->next == NULL){
            fail_msg("");
            return NULL;
        }
        child = child->next;
        child->parent = parent;
    }
    child->next = NULL;

    structure_node_add_area(parent, structure, classifier, result);
    
    return first_child;
}

static freesasa_structure_node *
structure_node_atom(const freesasa_structure *structure,
                    const freesasa_classifier *classifier,
                    const freesasa_result *result,
                    int atom_index)
{
    freesasa_structure_node *atom = 
        structure_node_new(freesasa_structure_atom_name(structure, atom_index));
                  
    if (atom == NULL) {
        fail_msg("");
        return NULL;
    } 
    atom->type = FREESASA_NODE_ATOM;
    atom->atom_prop = (struct atom_properties) {
        .is_polar = freesasa_classifier_class(classifier, freesasa_structure_atom_res_name(structure, atom_index), atom->name) == FREESASA_ATOM_POLAR,
        .is_bb = freesasa_atom_is_backbone(atom->name),
        .radius = freesasa_structure_atom_radius(structure, atom_index)
    };

    atom->area = malloc(sizeof(freesasa_nodearea));
    if (atom->area == NULL) {
        mem_fail();
        structure_node_free(atom);
        return NULL;
    }
    freesasa_atom_nodearea(atom->area, structure, result, classifier, atom_index);
    
    return atom;
}

static freesasa_structure_node *
structure_node_residue(const freesasa_structure *structure,
                       const freesasa_classifier *classifier,
                       const freesasa_result *result,
                       int residue_index)
{
    freesasa_structure_node *residue = NULL;
    const freesasa_nodearea *ref;
    int first, last;

    residue = structure_node_new(freesasa_structure_residue_name(structure, residue_index));
    
    if (residue == NULL) {
        fail_msg("");
        return NULL;
    }

    residue->type = FREESASA_NODE_RESIDUE;

    freesasa_structure_residue_atoms(structure, residue_index, &first, &last);
    residue->residue_prop.n_atoms = last - first + 1;

    residue->residue_prop.number = strdup(freesasa_structure_residue_number(structure, residue_index));
    if (residue->residue_prop.number == NULL) {
        mem_fail();
        goto cleanup;
    }

    ref = freesasa_classifier_residue_reference(classifier, residue->name);
    if (ref != NULL) {
        residue->residue_prop.reference = malloc(sizeof(freesasa_nodearea));
        if (residue->residue_prop.reference == NULL) {
            mem_fail();
            goto cleanup;
        }
        *residue->residue_prop.reference = *ref;
    } else {
        residue->residue_prop.reference = NULL;
    }

    if (structure_node_gen_children(residue, structure, classifier, result,
                                    first, last, structure_node_atom) == NULL) {
        goto cleanup;
    }

    return residue;

 cleanup:
    structure_node_free(residue);
    return NULL;
}

static freesasa_structure_node *
structure_node_chain(const freesasa_structure *structure,
                     const freesasa_classifier *classifier,
                     const freesasa_result *result,
                     int chain_index)
{
    const char *chains = freesasa_structure_chain_labels(structure);
    char name[2] = {chains[chain_index], '\0'};
    freesasa_structure_node *chain = NULL;
    int first_atom, last_atom, first_residue, last_residue;

    assert(strlen(chains) > chain_index);
    
    freesasa_structure_chain_atoms(structure, chains[chain_index], 
                                   &first_atom, &last_atom);

    chain = structure_node_new(name);
    if (chain == NULL) {
        fail_msg("");
        return NULL;
    }

    chain->type = FREESASA_NODE_CHAIN;
    freesasa_structure_chain_residues(structure, name[0],
                                      &first_residue, &last_residue);    
    chain->chain_prop.n_residues = last_residue - first_residue + 1;
    
    if (structure_node_gen_children(chain, structure, classifier, result,
                                    first_residue, last_residue,
                                    structure_node_residue) == NULL) {
        fail_msg("");
        structure_node_free(chain);
        return NULL;
    }

    return chain;
}

static freesasa_structure_node *
structure_node_structure(const freesasa_structure *structure,
                         const freesasa_classifier *classifier,
                         const freesasa_result *result,
                         int dummy_index)
{
    freesasa_structure_node *structure_node = NULL;

    structure_node = structure_node_new(freesasa_structure_chain_labels(structure));
    
    if (structure_node == NULL) {
        fail_msg("");
        return NULL;
    }

    structure_node->type = FREESASA_NODE_STRUCTURE;
    structure_node->structure_prop.n_chains = freesasa_structure_n_chains(structure);
    structure_node->structure_prop.chain_labels = strdup(freesasa_structure_chain_labels(structure));

    if (structure_node->structure_prop.chain_labels == NULL) {
        fail_msg("");
        goto cleanup;
    }
    
    if (structure_node_gen_children(structure_node, structure, classifier, result,
                                    0, freesasa_structure_n_chains(structure)-1,
                                    structure_node_chain) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return structure_node;
 cleanup:
    structure_node_free(structure_node);
    return NULL;
}

freesasa_structure_node *
freesasa_result2tree(const freesasa_result *result,
                     const freesasa_structure *structure,
                     const freesasa_classifier *classifier,
                     const char *name)
{
    if (classifier == NULL) classifier = &freesasa_default_classifier;
    freesasa_structure_node *root = structure_node_new(name);

    if (root == NULL) {
        goto cleanup;
    }

    root->type = FREESASA_NODE_ROOT;
    root->tree_prop.n_structures = 1;
    root->tree_prop.classified_by = strdup(classifier->name);
    if (root->tree_prop.classified_by == NULL) {
        goto cleanup;
    }
        
    if (structure_node_gen_children(root, structure, classifier, result,
                                    0, 0, structure_node_structure) == NULL) {
        goto cleanup;
    }

    return root;

 cleanup:
    structure_node_free(root);
    fail_msg("");
    return NULL;
}

int
freesasa_structure_node_join_trees(freesasa_structure_node *tree1,
                                   freesasa_structure_node **tree2)
{
    assert(tree1); assert(tree2); assert(*tree2);
    assert(tree1->type == FREESASA_NODE_ROOT);
    assert((*tree2)->type == FREESASA_NODE_ROOT);

    int ret = FREESASA_SUCCESS;

    if (strcmp(tree1->tree_prop.classified_by, (*tree2)->tree_prop.classified_by) != 0) {
        ret = freesasa_warn("Joining tree classified by different classifiers: '%s' and '%s', output may be inconsistent\n",
                            tree1->tree_prop.classified_by,
                            (*tree2)->tree_prop.classified_by);
    }

    freesasa_structure_node *child = tree1->children;
    if (child != NULL) {
        while (child->next) child = child->next;
    }
    child->next = (*tree2)->children;
    // tree1 takes over ownership, tree2 is invalidated.
    *tree2 = NULL;

    return ret;
}

int
freesasa_structure_node_free(freesasa_structure_node *root) 
{
    if (root) {
        if (root->parent)
            return fail_msg("Can't free node that isn't the root of its tree");
        structure_node_free(root);
    }
    return FREESASA_SUCCESS;
}

const freesasa_nodearea *
freesasa_structure_node_area(const freesasa_structure_node *node)
{
    assert(node->type != FREESASA_NODE_ROOT);
    return node->area;
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

const char*
freesasa_structure_node_classified_by(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_ROOT);
    return node->tree_prop.classified_by;
}

int
freesasa_structure_node_atom_is_polar(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->atom_prop.is_polar;
}

int
freesasa_structure_node_atom_is_mainchain(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->atom_prop.is_bb;
}

double
freesasa_structure_node_atom_radius(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->atom_prop.radius;
}

int
freesasa_structure_node_residue_n_atoms(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->residue_prop.n_atoms;
}

const char *
freesasa_structure_node_residue_number(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->residue_prop.number;
}

const freesasa_nodearea *
freesasa_structure_node_residue_reference(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->residue_prop.reference;
}

int
freesasa_structure_node_chain_n_residues(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_CHAIN);
    return node->chain_prop.n_residues;
}

int
freesasa_structure_node_structure_n_chains(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->structure_prop.n_chains;
}

const char *
freesasa_structure_node_structure_chain_labels(const freesasa_structure_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->structure_prop.chain_labels;
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
