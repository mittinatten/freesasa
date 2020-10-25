#if HAVE_CONFIG_H
#include <config.h>
#endif
#include <assert.h>
#include <errno.h>
#include <libxml/tree.h>
#include <libxml/xmlwriter.h>
#include <stdlib.h>
#include <string.h>

#include "freesasa_internal.h"
#include "pdb.h"

#ifndef FREESASA_XMLNS
#define FREESASA_XMLNS "freesasa"
#endif

static xmlNodePtr
nodearea2xml(const freesasa_nodearea *area,
             const char *name)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST name);
    char buf[20];

    sprintf(buf, "%.3f", area->total);
    if (xmlNewProp(xml_node, BAD_CAST "total", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%.3f", area->polar);
    if (xmlNewProp(xml_node, BAD_CAST "polar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%.3f", area->apolar);
    if (xmlNewProp(xml_node, BAD_CAST "apolar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%.3f", area->main_chain);
    if (xmlNewProp(xml_node, BAD_CAST "mainChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%.3f", area->side_chain);
    if (xmlNewProp(xml_node, BAD_CAST "sideChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return xml_node;

cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

static xmlNodePtr
atom2xml(const freesasa_node *node,
         int options)
{
    xmlNodePtr xml_node = NULL;
    const freesasa_nodearea *area = freesasa_node_area(node);
    const char *name = freesasa_node_name(node);
    char *trim_name, buf[20];

    assert(node);

    area = freesasa_node_area(node);
    name = freesasa_node_name(node);
    trim_name = malloc(strlen(name) + 1);

    if (!trim_name) {
        mem_fail();
        goto cleanup;
    }

    sscanf(name, "%s", trim_name);

    xml_node = xmlNewNode(NULL, BAD_CAST "atom");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST trim_name) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%.3f", area->total);
    if (xmlNewProp(xml_node, BAD_CAST "area", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%s", freesasa_node_atom_is_polar(node) == FREESASA_ATOM_POLAR ? "yes" : "no");
    if (xmlNewProp(xml_node, BAD_CAST "isPolar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%s", freesasa_atom_is_backbone(name) ? "yes" : "no");
    if (xmlNewProp(xml_node, BAD_CAST "isMainChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%.3f", freesasa_node_atom_radius(node));
    if (xmlNewProp(xml_node, BAD_CAST "radius", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    free(trim_name);
    return xml_node;

cleanup:
    free(trim_name);
    xmlFreeNode(xml_node);
    return NULL;
}

static xmlNodePtr
residue2xml(const freesasa_node *node,
            int options)
{
    xmlNodePtr xml_node = NULL, xml_area = NULL, xml_relarea = NULL;
    const char *name, *number;
    const freesasa_nodearea *abs, *reference;
    char *trim_number, *trim_name;
    freesasa_nodearea rel;

    assert(node);

    name = freesasa_node_name(node);
    number = freesasa_node_residue_number(node);
    abs = freesasa_node_area(node);
    reference = freesasa_node_residue_reference(node);
    trim_number = malloc(strlen(number) + 1);
    trim_name = malloc(strlen(name) + 1);

    if (!trim_number || !trim_name) {
        mem_fail();
        goto cleanup;
    }

    sscanf(number, "%s", trim_number);
    sscanf(name, "%s", trim_name);

    xml_node = xmlNewNode(NULL, BAD_CAST "residue");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST trim_name) == NULL ||
        xmlNewProp(xml_node, BAD_CAST "number", BAD_CAST trim_number) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = nodearea2xml(abs, "area");
    if (xml_area == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_node, xml_area) == NULL) {
        xmlFreeNode(xml_area);
        fail_msg("");
        goto cleanup;
    }

    if ((reference != NULL) && !(options & FREESASA_OUTPUT_SKIP_REL)) {
        freesasa_residue_rel_nodearea(&rel, abs, reference);
        xml_relarea = nodearea2xml(&rel, "relativeArea");
        if (xml_relarea == NULL) {
            fail_msg("");
            goto cleanup;
        }
        if (xmlAddChild(xml_node, xml_relarea) == NULL) {
            xmlFreeNode(xml_relarea);
            fail_msg("");
            goto cleanup;
        }
    }

    free(trim_name);
    free(trim_number);
    return xml_node;

cleanup:
    free(trim_name);
    free(trim_number);
    xmlFreeNode(xml_node);
    return NULL;
}

static xmlNodePtr
chain2xml(const freesasa_node *node,
          int options)
{
    xmlNodePtr xml_node = NULL, xml_area = NULL;
    char buf[20];

    xml_node = xmlNewNode(NULL, BAD_CAST "chain");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "label",
                   BAD_CAST freesasa_node_name(node)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%d", freesasa_node_chain_n_residues(node));
    if (xmlNewProp(xml_node, BAD_CAST "nResidues", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = nodearea2xml(freesasa_node_area(node), "area");
    if (xml_area == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_node, xml_area) == NULL) {
        fail_msg("");
        xmlFreeNode(xml_area);
        goto cleanup;
    }

    return xml_node;

cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

static xmlNodePtr
selection2xml(const freesasa_selection *selection)
{
    xmlNodePtr xml_selection = xmlNewNode(NULL, BAD_CAST "selection");
    char buf[20];

    sprintf(buf, "%.3f", freesasa_selection_area(selection));

    if (xml_selection == NULL) {
        fail_msg("");
    } else {
        if (xmlNewProp(xml_selection, BAD_CAST "name", BAD_CAST freesasa_selection_name(selection)) == NULL ||
            xmlNewProp(xml_selection, BAD_CAST "area", BAD_CAST buf) == NULL) {
            fail_msg("");
            xmlFreeNode(xml_selection);
            xml_selection = NULL;
        }
    }

    return xml_selection;
}

static xmlNodePtr
add_selections_to_xmlNode(const freesasa_selection **selections,
                          xmlNodePtr parent)
{
    xmlNodePtr xml_selection;

    while (*selections) {
        xml_selection = selection2xml(*selections);
        if (xml_selection == NULL) {
            fail_msg("");
            return NULL;
        }

        if (xmlAddChild(parent, xml_selection) == NULL) {
            fail_msg("");
            xmlFreeNode(xml_selection);
            return NULL;
        }
        ++selections;
    }
    return parent;
}

static xmlNodePtr
structure2xml(const freesasa_node *node,
              int options)
{
    xmlNodePtr xml_node = NULL, xml_area = NULL;
    char buf[20];
    const freesasa_selection **selections;

    assert(node);

    selections = freesasa_node_structure_selections(node);

    xml_node = xmlNewNode(NULL, BAD_CAST "structure");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "chains", BAD_CAST freesasa_node_structure_chain_labels(node)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%d", freesasa_node_structure_model(node));
    if (xmlNewProp(xml_node, BAD_CAST "model", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = nodearea2xml(freesasa_node_area(node), "area");
    if (xml_area == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_node, xml_area) == NULL) {
        fail_msg("");
        xmlFreeNode(xml_area);
        goto cleanup;
    }

    if (selections) {
        if (add_selections_to_xmlNode(selections, xml_node) == NULL) {
            fail_msg("");
            goto cleanup;
        }
    }

    return xml_node;

cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

/* Root is not converted to xmlNodePtr, we skip immediately to children */
static int
node2xml(xmlNodePtr *xml_node,
         freesasa_node *node,
         int exclude_type,
         int options)
{
    freesasa_node *child;
    xmlNodePtr xml_child = NULL;

    assert(xml_node);
    assert(node);

    child = freesasa_node_children(node);
    *xml_node = NULL;

    if (freesasa_node_type(node) == exclude_type) return FREESASA_SUCCESS;

    switch (freesasa_node_type(node)) {
    case FREESASA_NODE_STRUCTURE:
        *xml_node = structure2xml(node, options);
        break;
    case FREESASA_NODE_CHAIN:
        *xml_node = chain2xml(node, options);
        break;
    case FREESASA_NODE_RESIDUE:
        *xml_node = residue2xml(node, options);
        break;
    case FREESASA_NODE_ATOM:
        *xml_node = atom2xml(node, options);
        break;
    case FREESASA_NODE_ROOT:
    default:
        assert(0 && "tree illegal");
    }
    if (*xml_node == NULL)
        return fail_msg("error creating XML-node");

    /* simplify? */
    while (child != NULL) {
        if (node2xml(&xml_child, child, exclude_type, options) == FREESASA_FAIL) {
            return fail_msg("");
        }

        if (xml_child != NULL &&
            xmlAddChild(*xml_node, xml_child) == NULL) {
            xmlFreeNode(*xml_node);
            xmlFreeNode(xml_child);
            return fail_msg("");
        }
        child = freesasa_node_next(child);
    }

    return FREESASA_SUCCESS;
}

static xmlNodePtr
parameters2xml(const freesasa_parameters *p)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "parameters");
    char buf[20];

    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "algorithm", BAD_CAST freesasa_alg_name(p->alg)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", p->probe_radius);
    if (xmlNewProp(xml_node, BAD_CAST "probeRadius", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    switch (p->alg) {
    case FREESASA_SHRAKE_RUPLEY:
        sprintf(buf, "%d", p->shrake_rupley_n_points);
        break;
    case FREESASA_LEE_RICHARDS:
        sprintf(buf, "%d", p->lee_richards_n_slices);
        break;
    default:
        assert(0);
        break;
    }
    if (xmlNewProp(xml_node, BAD_CAST "resolution", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return xml_node;

cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

static xmlNodePtr
xml_result(freesasa_node *result,
           int options)
{
    xmlNodePtr xml_result_node = NULL, xml_structure = NULL, xml_param = NULL;
    freesasa_node *child = NULL;
    const freesasa_parameters *parameters;
    int exclude_type = FREESASA_NODE_NONE;

    assert(freesasa_node_type(result) == FREESASA_NODE_RESULT);

    parameters = freesasa_node_result_parameters(result);

    if (options & FREESASA_OUTPUT_STRUCTURE) exclude_type = FREESASA_NODE_CHAIN;
    if (options & FREESASA_OUTPUT_CHAIN) exclude_type = FREESASA_NODE_RESIDUE;
    if (options & FREESASA_OUTPUT_RESIDUE) exclude_type = FREESASA_NODE_ATOM;

    xml_result_node = xmlNewNode(NULL, BAD_CAST "result");
    if (xml_result_node == NULL) {
        fail_msg("");
        return NULL;
    }

    xml_param = parameters2xml(parameters);
    if (xml_param == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_result_node, xml_param) == NULL) {
        fail_msg("");
        xmlFreeNode(xml_param);
        goto cleanup;
    }

    if (xmlNewProp(xml_result_node, BAD_CAST "classifier",
                   BAD_CAST freesasa_node_classified_by(result)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlNewProp(xml_result_node, BAD_CAST "input",
                   BAD_CAST freesasa_node_name(result)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    child = freesasa_node_children(result);
    assert(child);

    while (child) {
        if (node2xml(&xml_structure, child, exclude_type, options) == FREESASA_FAIL) {
            fail_msg("");
            goto cleanup;
        }
        if (xmlAddChild(xml_result_node, xml_structure) == NULL) {
            fail_msg("");
            goto cleanup;
        }
        child = freesasa_node_next(child);
    };

    return xml_result_node;
cleanup:
    xmlFreeNode(xml_structure);
    xmlFreeNode(xml_result_node);
    return NULL;
}

int freesasa_write_xml(FILE *output,
                       freesasa_node *root,
                       int options)
{
    freesasa_node *child = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr xml_root = NULL, xml_result_node = NULL;
    xmlNsPtr ns = NULL;
    xmlBufferPtr buf = NULL;
    xmlTextWriterPtr writer = NULL;
    int ret = FREESASA_FAIL;

    assert(freesasa_node_type(root) == FREESASA_NODE_ROOT);

    doc = xmlNewDoc(BAD_CAST "1.0");
    if (doc == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_root = xmlNewNode(NULL, BAD_CAST "results");
    if (xml_root == NULL) {
        fail_msg("");
        goto cleanup;
    }

    ns = xmlNewNs(xml_root, BAD_CAST FREESASA_XMLNS, NULL);
    if (ns == NULL) {
        fail_msg("");
        xmlFreeNode(xml_root);
        goto cleanup;
    }

    xmlDocSetRootElement(doc, xml_root);

    /* global attributes */
    if (xmlNewProp(xml_root, BAD_CAST "source", BAD_CAST freesasa_string) == NULL) {
        fail_msg("");
        goto cleanup;
    }
    if (xmlNewProp(xml_root, BAD_CAST "lengthUnit", BAD_CAST "Ångström") == NULL) {
        fail_msg("");
        goto cleanup;
    }

    child = freesasa_node_children(root);
    while (child) {
        xml_result_node = xml_result(child, options);
        if (xml_result_node == NULL) {
            fail_msg("");
            goto cleanup;
        }
        if (xmlAddChild(xml_root, xml_result_node) == NULL) {
            fail_msg("");
            xmlFreeNode(xml_result_node);
            goto cleanup;
        }
        child = freesasa_node_next(child);
    }

    buf = xmlBufferCreate();
    if (buf == NULL) {
        fail_msg("");
        goto cleanup;
    }

    writer = xmlNewTextWriterMemory(buf, 0);
    if (writer == NULL) {
        xmlBufferFree(buf);
        fail_msg("");
        goto cleanup;
    }

    if (xmlTextWriterStartDocument(writer, XML_DEFAULT_VERSION,
                                   xmlGetCharEncodingName(XML_CHAR_ENCODING_UTF8), NULL) == -1) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlTextWriterFlush(writer) == -1) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlNodeDump(buf, doc, xml_root, 0, 1) == 0) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlTextWriterEndDocument(writer) == -1) {
        fail_msg("");
        goto cleanup;
    }

    fprintf(output, "%s", (const char *)buf->content);
    fflush(output);
    if (ferror(output)) {
        fail_msg(strerror(errno));
        goto cleanup;
    }

    ret = FREESASA_SUCCESS;

cleanup:
    xmlFreeDoc(doc);
    xmlFreeTextWriter(writer);
    return ret;
}
