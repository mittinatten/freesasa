#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "classifier.h"
#include "freesasa_internal.h"

struct atom_properties {
    int is_polar;
    int is_bb;
    double radius;
    char *pdb_line;
    char *chain;
    char *res_number;
    char *res_name;
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
    size_t cif_ref;
    freesasa_result *result;
    freesasa_selection **selection; // NULL terminated array
};

struct result_properties {
    char *classified_by;
    freesasa_parameters parameters;
    int n_structures;
};

struct freesasa_node {
    char *name;
    freesasa_nodetype type;
    union properties {
        struct atom_properties atom;
        struct residue_properties residue;
        struct chain_properties chain;
        struct structure_properties structure;
        struct result_properties result;
    } properties;
    freesasa_nodearea *area;
    freesasa_node *parent;
    freesasa_node *children;
    freesasa_node *next;
};

const freesasa_nodearea freesasa_nodearea_null = {NULL, 0, 0, 0, 0, 0, 0};

static freesasa_node *
node_new(const char *name)
{
    freesasa_node *node = malloc(sizeof(freesasa_node));

    if (node == NULL) {
        goto memerr;
    }

    node->name = NULL;
    node->type = FREESASA_NODE_ATOM;
    node->area = NULL;
    node->parent = NULL;
    node->children = NULL;
    node->next = NULL;

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
node_free(freesasa_node *node)
{
    freesasa_node *current = NULL, *next = NULL;
    freesasa_selection **sel = NULL;

    if (node != NULL) {
        current = node->children;
        while (current) {
            next = current->next;
            node_free(current);
            current = next;
        }
        free(node->name);
        free(node->area);

        switch (node->type) {
        case FREESASA_NODE_ATOM:
            free(node->properties.atom.pdb_line);
            free(node->properties.atom.res_name);
            free(node->properties.atom.res_number);
            free(node->properties.atom.chain);
            break;
        case FREESASA_NODE_RESIDUE:
            free(node->properties.residue.reference);
            free(node->properties.residue.number);
            break;
        case FREESASA_NODE_STRUCTURE:
            free(node->properties.structure.chain_labels);
            freesasa_result_free(node->properties.structure.result);
            sel = node->properties.structure.selection;
            if (sel) {
                while (*sel) {
                    freesasa_selection_free(*sel);
                    ++sel;
                }
            }
            free(node->properties.structure.selection);
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

typedef freesasa_node *(*node_generator)(const freesasa_structure *,
                                         const freesasa_result *,
                                         int index);

static int
node_add_area(freesasa_node *node,
              const freesasa_structure *structure,
              const freesasa_result *result)
{
    freesasa_node *child = NULL;

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

static freesasa_node *
node_gen_children(freesasa_node *parent,
                  const freesasa_structure *structure,
                  const freesasa_result *result,
                  int first,
                  int last,
                  node_generator ng)
{
    int i;
    freesasa_node *child, *first_child;

    first_child = ng(structure, result, first);

    if (first_child == NULL) {
        fail_msg("");
        return NULL;
    }

    first_child->parent = parent;
    child = parent->children = first_child;

    for (i = first + 1; i <= last; ++i) {
        child->next = ng(structure, result, i);
        if (child->next == NULL) {
            fail_msg("");
            return NULL;
        }
        child = child->next;
        child->parent = parent;
    }
    child->next = NULL;

    node_add_area(parent, structure, result);

    return first_child;
}

static freesasa_node *
node_atom(const freesasa_structure *structure,
          const freesasa_result *result,
          int atom_index)
{
    freesasa_node *atom =
        node_new(freesasa_structure_atom_name(structure, atom_index));
    const char *line;

    if (atom == NULL) {
        fail_msg("");
        return NULL;
    }

    atom->type = FREESASA_NODE_ATOM;
    atom->properties.atom.pdb_line = NULL;
    atom->properties.atom.res_number = NULL;
    atom->properties.atom.res_name = NULL;
    atom->properties.atom.is_polar = freesasa_structure_atom_class(structure, atom_index) == FREESASA_ATOM_POLAR;
    atom->properties.atom.is_bb = freesasa_atom_is_backbone(atom->name);
    atom->properties.atom.radius = freesasa_structure_atom_radius(structure, atom_index);

    atom->properties.atom.chain = strdup(freesasa_structure_atom_chain_lcl(structure, atom_index));
    if (atom->properties.atom.chain == NULL) {
        mem_fail();
        goto cleanup;
    }

    atom->properties.atom.res_number = strdup(freesasa_structure_atom_res_number(structure, atom_index));
    if (atom->properties.atom.res_number == NULL) {
        mem_fail();
        goto cleanup;
    }

    atom->properties.atom.res_name = strdup(freesasa_structure_atom_res_name(structure, atom_index));
    if (atom->properties.atom.res_name == NULL) {
        mem_fail();
        goto cleanup;
    }

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
    node_free(atom);
    return NULL;
}

static freesasa_node *
node_residue(const freesasa_structure *structure,
             const freesasa_result *result,
             int residue_index)
{
    freesasa_node *residue = NULL;
    const freesasa_nodearea *ref;
    int first, last;

    residue = node_new(freesasa_structure_residue_name(structure, residue_index));

    if (residue == NULL) {
        fail_msg("");
        return NULL;
    }

    residue->type = FREESASA_NODE_RESIDUE;

    freesasa_structure_residue_atoms(structure, residue_index, &first, &last);
    residue->properties.residue.n_atoms = last - first + 1;
    residue->properties.residue.reference = NULL;

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

    if (node_gen_children(residue, structure, result, first,
                          last, node_atom) == NULL) {
        goto cleanup;
    }

    return residue;

cleanup:
    node_free(residue);
    return NULL;
}

static freesasa_node *
node_chain(const freesasa_structure *structure,
           const freesasa_result *result,
           int chain_index)
{
    const char *name = freesasa_structure_chain_label(structure, chain_index);
    freesasa_node *chain = NULL;
    int first_atom, last_atom, first_residue, last_residue;

    freesasa_structure_chain_atoms_lcl(structure, name,
                                       &first_atom, &last_atom);

    chain = node_new(name);
    if (chain == NULL) {
        fail_msg("");
        return NULL;
    }

    chain->type = FREESASA_NODE_CHAIN;
    freesasa_structure_chain_residues_lcl(structure, name,
                                          &first_residue, &last_residue);
    chain->properties.chain.n_residues = last_residue - first_residue + 1;

    if (node_gen_children(chain, structure, result,
                          first_residue, last_residue,
                          node_residue) == NULL) {
        fail_msg("");
        node_free(chain);
        return NULL;
    }

    return chain;
}

static freesasa_node *
node_structure(const freesasa_structure *structure,
               const freesasa_result *result,
               int dummy_index)
{
    freesasa_node *node = NULL;
    node = node_new(freesasa_structure_chain_labels(structure));

    if (node == NULL) {
        fail_msg("");
        return NULL;
    }

    node->type = FREESASA_NODE_STRUCTURE;
    node->properties.structure.n_chains = freesasa_structure_n_chains(structure);
    node->properties.structure.n_atoms = freesasa_structure_n(structure);
    node->properties.structure.result = NULL;
    node->properties.structure.selection = NULL;
    node->properties.structure.chain_labels = strdup(freesasa_structure_chain_labels(structure));
    node->properties.structure.model = freesasa_structure_model(structure);
    node->properties.structure.cif_ref = freesasa_structure_cif_ref(structure);

    if (node->properties.structure.chain_labels == NULL) {
        mem_fail();
        goto cleanup;
    }

    node->properties.structure.result = freesasa_result_clone(result);

    if (node->properties.structure.result == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (node_gen_children(node, structure, result, 0,
                          freesasa_structure_n_chains(structure) - 1,
                          node_chain) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return node;
cleanup:
    node_free(node);
    return NULL;
}

freesasa_node *
freesasa_tree_new(void)
{
    freesasa_node *tree = node_new(NULL);
    if (tree != NULL) {
        tree->type = FREESASA_NODE_ROOT;
    }
    return tree;
}

freesasa_node *
freesasa_tree_init(const freesasa_result *result,
                   const freesasa_structure *structure,
                   const char *name)
{
    freesasa_node *tree = node_new(NULL);

    tree->type = FREESASA_NODE_ROOT;

    if (tree == NULL) {
        fail_msg("");
    } else if (freesasa_tree_add_result(tree, result, structure, name) == FREESASA_FAIL) {
        fail_msg("");
        freesasa_node_free(tree);
        tree = NULL;
    }

    return tree;
}

int freesasa_tree_add_result(freesasa_node *tree,
                             const freesasa_result *result,
                             const freesasa_structure *structure,
                             const char *name)
{
    freesasa_node *node = node_new(name);

    if (node == NULL) {
        goto cleanup;
    }

    node->type = FREESASA_NODE_RESULT;
    node->properties.result.n_structures = 1;
    node->properties.result.parameters = result->parameters;
    node->properties.result.classified_by = strdup(freesasa_structure_classifier_name(structure));

    if (node->properties.result.classified_by == NULL) {
        mem_fail();
        goto cleanup;
    }

    if (node_gen_children(node, structure, result, 0, 0,
                          node_structure) == NULL) {
        goto cleanup;
    }

    node->next = tree->children;
    tree->children = node;

    return FREESASA_SUCCESS;

cleanup:
    node_free(node);
    fail_msg("");
    return FREESASA_FAIL;
}

int freesasa_tree_join(freesasa_node *tree1,
                       freesasa_node **tree2)
{
    freesasa_node *child;

    assert(tree1);
    assert(tree2);
    assert(*tree2);
    assert(tree1->type == FREESASA_NODE_ROOT);
    assert((*tree2)->type == FREESASA_NODE_ROOT);

    child = tree1->children;

    if (child != NULL) {
        while (child->next)
            child = child->next;
        child->next = (*tree2)->children;
    } else {
        tree1->children = (*tree2)->children;
    }
    // tree1 takes over ownership, tree2 is invalidated.
    free(*tree2);
    *tree2 = NULL;

    return FREESASA_SUCCESS;
}

int freesasa_node_free(freesasa_node *root)
{
    if (root) {
        if (root->parent)
            return fail_msg("can't free node that isn't the root of its tree");
        node_free(root);
    }
    return FREESASA_SUCCESS;
}

const freesasa_nodearea *
freesasa_node_area(const freesasa_node *node)
{
    assert(node->type != FREESASA_NODE_ROOT);
    return node->area;
}

freesasa_node *
freesasa_node_children(freesasa_node *node)
{
    return node->children;
}

freesasa_node *
freesasa_node_next(freesasa_node *node)
{
    return node->next;
}

freesasa_node *
freesasa_node_parent(freesasa_node *node)
{
    return node->parent;
}

freesasa_nodetype
freesasa_node_type(const freesasa_node *node)
{
    return node->type;
}

const char *
freesasa_node_name(const freesasa_node *node)
{
    return node->name;
}

const char *
freesasa_node_classified_by(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_RESULT);
    return node->properties.result.classified_by;
}

int freesasa_node_atom_is_polar(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.is_polar;
}

int freesasa_node_atom_is_mainchain(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.is_bb;
}

double
freesasa_node_atom_radius(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.radius;
}

const char *
freesasa_node_atom_pdb_line(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.pdb_line;
}

const char *
freesasa_node_atom_residue_number(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.res_number;
}

const char *
freesasa_node_atom_residue_name(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.res_name;
}

char *freesasa_node_atom_chain(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_ATOM);
    return node->properties.atom.chain;
}

int freesasa_node_residue_n_atoms(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->properties.residue.n_atoms;
}

const char *
freesasa_node_residue_number(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->properties.residue.number;
}

const freesasa_nodearea *
freesasa_node_residue_reference(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_RESIDUE);
    return node->properties.residue.reference;
}

int freesasa_node_chain_n_residues(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_CHAIN);
    return node->properties.chain.n_residues;
}

int freesasa_node_structure_n_chains(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.n_chains;
}

int freesasa_node_structure_n_atoms(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.n_atoms;
}

int freesasa_node_structure_model(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.model;
}

const char *
freesasa_node_structure_chain_labels(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.chain_labels;
}

const freesasa_result *
freesasa_node_structure_result(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.result;
}

size_t
freesasa_node_structure_cif_ref(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return node->properties.structure.cif_ref;
}

int freesasa_node_structure_add_selection(freesasa_node *node,
                                          const freesasa_selection *selection)
{
    int n_selections = 0;
    freesasa_selection **sel, **sel2;

    assert(node->type == FREESASA_NODE_STRUCTURE);

    sel = node->properties.structure.selection;

    // count number of selections
    if (sel) {
        sel2 = sel;
        while (*sel2) {
            ++n_selections;
            ++sel2;
        }
    }

    sel = realloc(sel, sizeof(freesasa_selection *) * (n_selections + 2));
    if (sel == NULL) {
        return mem_fail();
    }

    sel[n_selections] = freesasa_selection_clone(selection);
    if (sel[n_selections] == NULL) {
        return fail_msg("");
    }
    sel[n_selections + 1] = NULL;
    node->properties.structure.selection = sel;

    return FREESASA_SUCCESS;
}

const freesasa_selection **
freesasa_node_structure_selections(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_STRUCTURE);
    return (const freesasa_selection **)node->properties.structure.selection;
}

const freesasa_parameters *
freesasa_node_result_parameters(const freesasa_node *node)
{
    assert(node->type == FREESASA_NODE_RESULT);
    return &node->properties.result.parameters;
}

int freesasa_atom_nodearea(freesasa_nodearea *area,
                           const freesasa_structure *structure,
                           const freesasa_result *result,
                           int atom_index)
{
    double a = result->sasa[atom_index];
    *area = freesasa_nodearea_null;

    area->total = a;

    if (freesasa_atom_is_backbone(freesasa_structure_atom_name(structure, atom_index)))
        area->main_chain = a;
    else
        area->side_chain = a;

    switch (freesasa_structure_atom_class(structure, atom_index)) {
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

void freesasa_add_nodearea(freesasa_nodearea *sum,
                           const freesasa_nodearea *term)
{
    sum->total += term->total;
    sum->side_chain += term->side_chain;
    sum->main_chain += term->main_chain;
    sum->polar += term->polar;
    sum->apolar += term->apolar;
    sum->unknown += term->unknown;
}

void freesasa_range_nodearea(freesasa_nodearea *area,
                             const freesasa_structure *structure,
                             const freesasa_result *result,
                             int first_atom,
                             int last_atom)
{
    int i;
    freesasa_nodearea term = freesasa_nodearea_null;

    assert(area);
    assert(structure);
    assert(result);
    assert(first_atom <= last_atom);

    for (i = first_atom; i <= last_atom; ++i) {
        freesasa_atom_nodearea(&term, structure, result, i);
        freesasa_add_nodearea(area, &term);
    }
}
