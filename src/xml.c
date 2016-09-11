#if HAVE_CONFIG_H
  #include <config.h>
#endif
#include <libxml/tree.h>
#include <libxml/xmlwriter.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include "freesasa_internal.h"

static xmlNodePtr
xml_nodearea(const freesasa_nodearea *area, const char *name)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST name);
    char buf[20];

    sprintf(buf, "%f", area->total);
    if (xmlNewProp(xml_node, BAD_CAST "total", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", area->polar);
    if (xmlNewProp(xml_node, BAD_CAST "polar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }
    
    sprintf(buf, "%f", area->apolar);
    if (xmlNewProp(xml_node, BAD_CAST "apolar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", area->main_chain);
    if (xmlNewProp(xml_node, BAD_CAST "mainChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", area->side_chain);
    if (xmlNewProp(xml_node, BAD_CAST "sideChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return xml_node;

 cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

xmlNodePtr
freesasa_xml_atom(const freesasa_result_node *node, int options)
{
    assert(node);
    xmlNodePtr xml_node = NULL;
    const freesasa_nodearea *area = freesasa_result_node_area(node);
    const char *name = freesasa_result_node_name(node);
    int n_len = strlen(name);
    char trim_name[n_len+1], buf[20];

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

    sprintf(buf, "%f", area->total);
    if (xmlNewProp(xml_node, BAD_CAST "area", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%s", freesasa_result_node_atom_is_polar(node) == FREESASA_ATOM_POLAR ? "yes" : "no");
    if (xmlNewProp(xml_node, BAD_CAST "isPolar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%s", freesasa_atom_is_backbone(name) ? "yes" : "no");
    if (xmlNewProp(xml_node, BAD_CAST "isMainChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", freesasa_result_node_atom_radius(node));
    if (xmlNewProp(xml_node, BAD_CAST "radius", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return xml_node;

 cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

xmlNodePtr
freesasa_xml_residue(const freesasa_result_node *node, int options)
{
    assert(node);
    xmlNodePtr xml_node = NULL, xml_area = NULL, xml_relarea = NULL;
    const char *name = freesasa_result_node_name(node), *number;
    const freesasa_nodearea *abs = freesasa_result_node_area(node),
        *reference = freesasa_result_node_residue_reference(node);
    freesasa_nodearea rel;

    number = freesasa_result_node_residue_number(node);

    int n_len = strlen(number);
    char trim_number[n_len+1];
    sscanf(number, "%s", trim_number);
    xml_node = xmlNewNode(NULL, BAD_CAST "residue");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST name) == NULL ||
        xmlNewProp(xml_node, BAD_CAST "number", BAD_CAST trim_number) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = xml_nodearea(abs, "area");
    if (xml_area == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_node, xml_area) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if ((reference != NULL) && !(options & FREESASA_OUTPUT_SKIP_REL)) {
        freesasa_residue_rel_nodearea(&rel, abs, reference);
        xml_relarea = xml_nodearea(&rel, "relativeArea");
        if (xml_relarea == NULL) {
            fail_msg("");
            goto cleanup;
        }
        if (xmlAddChild(xml_node, xml_relarea) == NULL) {
            fail_msg("");
            goto cleanup;
        }
    }

    return xml_node;

 cleanup:
    xmlFreeNode(xml_node);
    return NULL;
}

xmlNodePtr
freesasa_xml_chain(const freesasa_result_node *node, int options)
{
    xmlNodePtr xml_node = NULL, xml_area = NULL;
    const char *name = freesasa_result_node_name(node);
    char buf[20];

    xml_node = xmlNewNode(NULL, BAD_CAST "chain");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "label",
                   BAD_CAST freesasa_result_node_name(node)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%d", freesasa_result_node_chain_n_residues(node));
    if (xmlNewProp(xml_node, BAD_CAST "nResidues", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = xml_nodearea(freesasa_result_node_area(node), "area");
    if (xml_area == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_node, xml_area) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return xml_node;

 cleanup:
    xmlFreeNode(xml_area);
    xmlFreeNode(xml_node);
    return NULL;
}

xmlNodePtr
freesasa_xml_structure(const freesasa_result_node *node, int options)
{
    assert(node);
    xmlNodePtr xml_node = NULL, xml_area = NULL;

    xml_node = xmlNewNode(NULL, BAD_CAST "structure");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "chains", BAD_CAST freesasa_result_node_structure_chain_labels(node)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = xml_nodearea(freesasa_result_node_area(node), "area");
    if (xmlAddChild(xml_node, xml_area) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    return xml_node;

 cleanup:
    xmlFreeNode(xml_node);
    xmlFreeNode(xml_area);
    return NULL;
}

// Root is not converted to xmlNodePtr, we skip immediately to children
int
freesasa_node2xml(xmlNodePtr *xml_node,
                  const freesasa_result_node *node,
                  int exclude_type,
                  int options)
{
    assert(xml_node);
    assert(node);
    const freesasa_result_node *child = freesasa_result_node_children(node);
    xmlNodePtr xml_child = NULL;
    *xml_node = NULL;

    if (freesasa_result_node_type(node) == exclude_type) return FREESASA_SUCCESS;

    switch (freesasa_result_node_type(node)) {
    case FREESASA_NODE_STRUCTURE:
        *xml_node = freesasa_xml_structure(node, options);
        break;
    case FREESASA_NODE_CHAIN:
        *xml_node = freesasa_xml_chain(node, options);
        break;
    case FREESASA_NODE_RESIDUE:
        *xml_node = freesasa_xml_residue(node, options);
        break;
    case FREESASA_NODE_ATOM:
        *xml_node = freesasa_xml_atom(node, options);
        break;
    case FREESASA_NODE_ROOT:
    default:
        assert(0 && "Tree illegal");
    }
    if (*xml_node == NULL)
        return fail_msg("");

    // simplify?
    while (child != NULL) {
        if (freesasa_node2xml(&xml_child, child, exclude_type, options) == FREESASA_FAIL)
            return fail_msg("");

        if (xml_child != NULL &&
            xmlAddChild(*xml_node, xml_child) == NULL) {
            xmlFreeNode(xml_child);
            return fail_msg("");
        }
        child = freesasa_result_node_next(child);
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

    switch(p->alg) {
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
xml_result(const freesasa_result_node *result,
           const freesasa_parameters *parameters,
           int options)
{
    assert(freesasa_result_node_type(result) == FREESASA_NODE_RESULT);
    xmlNodePtr xml_result_node = NULL, xml_structure = NULL, xml_param = NULL;
    const freesasa_result_node *child = NULL;
    int exclude_type = FREESASA_NODE_NONE;

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
                   BAD_CAST freesasa_result_node_classified_by(result)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlNewProp(xml_result_node, BAD_CAST "input",
                   BAD_CAST freesasa_result_node_name(result)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    child = freesasa_result_node_children(result);
    assert(child);

    while(child) {
        if (freesasa_node2xml(&xml_structure, child, exclude_type, options) == FREESASA_FAIL) {
            fail_msg("");
            goto cleanup;
        }
        if (xmlAddChild(xml_result_node, xml_structure) == NULL) {
            fail_msg("");
            goto cleanup;
        }
        child = freesasa_result_node_next(child);
    };

    return xml_result_node;
 cleanup:
    xmlFreeNode(xml_result_node);
    return NULL;
}

int
freesasa_write_xml(FILE *output,
                   const freesasa_result_node *root,
                   const freesasa_parameters *parameters,
                   int options)
{
    assert(freesasa_result_node_type(root) == FREESASA_NODE_ROOT);

    const freesasa_result_node *child = NULL;
    xmlDocPtr doc = NULL;
    xmlNodePtr xml_root = NULL, xml_result_node = NULL;
    xmlNsPtr ns = NULL;
    xmlBufferPtr buf = NULL;
    xmlTextWriterPtr writer = NULL;
    int ret = FREESASA_FAIL;

    if (parameters == NULL) parameters = &freesasa_default_parameters;

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

    ns = xmlNewNs(xml_root, BAD_CAST "http://freesasa.github.io/", NULL);
    buf = xmlBufferCreate();
    if (ns == NULL || buf == NULL) {
        fail_msg("");
        // this will later be stored with the doc, so can't be freed twice in cleanup
        xmlFreeNs(ns);
        goto cleanup;
    }

    xmlDocSetRootElement(doc, xml_root);

    // global attributes
    if (xmlNewProp(xml_root, BAD_CAST "source", BAD_CAST freesasa_string) == NULL) {
        fail_msg("");
        goto cleanup;
    }
    if (xmlNewProp(xml_root, BAD_CAST "lengthUnit", BAD_CAST "Ångström") == NULL) {
        fail_msg("");
        goto cleanup;
    }

    child = freesasa_result_node_children(root);
    while (child) {
        xml_result_node = xml_result(child, parameters, options);
        if (xml_result_node == NULL) {
            fail_msg("");
            goto cleanup;
        }
        if (xmlAddChild(xml_root, xml_result_node) == NULL) {
            fail_msg("");
            xmlFreeNode(xml_result_node);
            goto cleanup;
        }
        child = freesasa_result_node_next(child);
    }

    writer = xmlNewTextWriterMemory(buf, 0);
    if (writer == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlTextWriterStartDocument(writer, XML_DEFAULT_VERSION,
                                   xmlGetCharEncodingName(XML_CHAR_ENCODING_UTF8), NULL)
        == -1) {
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
    
    fprintf(output, "%s", (const char*) buf->content);
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
