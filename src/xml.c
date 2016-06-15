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
xml_subarea(const freesasa_subarea *area, const char *name)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, name);
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
freesasa_xml_atom(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = NULL;
    const freesasa_subarea *area = freesasa_structure_node_area(node);
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    int first, last;
    double radius;
    const char *name = freesasa_structure_node_name(node);
    int n_len = strlen(name);
    char trim_name[n_len+1], buf[20];
    const char *resn = freesasa_structure_node_name(freesasa_structure_node_parent(node));
    int is_polar;
    int is_bb = freesasa_atom_is_backbone(name);
    double sasa = area->total;

    sscanf(name, "%s", trim_name);
    freesasa_structure_node_atoms(node, &first, &last);
    radius = freesasa_structure_atom_radius(structure, first);
    is_polar = freesasa_structure_atom_class(structure, first);

    xml_node = xmlNewNode(NULL, BAD_CAST "atom");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST trim_name) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", sasa);
    if (xmlNewProp(xml_node, BAD_CAST "area", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%s", is_polar ? "yes" : "no");
    if (xmlNewProp(xml_node, BAD_CAST "isPolar", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%s", is_bb ? "yes" : "no");
    if (xmlNewProp(xml_node, BAD_CAST "isMainChain", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%f", radius);
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
freesasa_xml_residue(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = NULL, xml_area = NULL, xml_relarea = NULL;
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node), *number;
    const freesasa_subarea *abs = freesasa_structure_node_area(node),
        *reference = freesasa_structure_node_residue_reference(node);
    freesasa_subarea rel;
    int first, last;

    freesasa_structure_node_atoms(node, &first, &last);
    number = freesasa_structure_atom_res_number(structure, first);

    int n_len = strlen(number);
    char trim_number[n_len+1];
    sscanf(number, "%s", trim_number);
    xml_node = xmlNewNode(NULL, "residue");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST name) == NULL ||
        xmlNewProp(xml_node, BAD_CAST "number", BAD_CAST trim_number) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = xml_subarea(abs, "area");
    if (xml_area == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_node, xml_area) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (reference != NULL) {
        freesasa_residue_rel_subarea(&rel, abs, reference);
        xml_relarea = xml_subarea(&rel, "relativeArea");
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
freesasa_xml_chain(const freesasa_structure_node *node)
{
    xmlNodePtr xml_node = NULL, xml_area = NULL;
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node);
    char buf[20];
    int first, last;
    freesasa_structure_chain_residues(structure, name[0], &first, &last);

    xml_node = xmlNewNode(NULL, BAD_CAST "chain");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    if (xmlNewProp(xml_node, BAD_CAST "label",
                   BAD_CAST freesasa_structure_node_name(node)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    sprintf(buf, "%d", last - first + 1);
    if (xmlNewProp(xml_node, BAD_CAST "nResidues", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = xml_subarea(freesasa_structure_node_area(node), "area");
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
freesasa_xml_structure(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = NULL, xml_area = NULL;
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    char buf[20];

    xml_node = xmlNewNode(NULL, BAD_CAST "structure");
    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

    sprintf(buf, "%d", freesasa_structure_n_chains(structure));
    if (xmlNewProp(xml_node, BAD_CAST "nChains", BAD_CAST buf) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_area = xml_subarea(freesasa_structure_node_area(node), "area");
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

int
freesasa_node2xml(xmlNodePtr *xml_node, const freesasa_structure_node *node, int exclude_type)
{
    assert(xml_node);
    assert(node);
    const freesasa_structure_node *child = freesasa_structure_node_children(node);
    xmlNodePtr xml_child = NULL;
    *xml_node == NULL;

    if (freesasa_structure_node_type(node) == exclude_type) return FREESASA_SUCCESS;

    switch (freesasa_structure_node_type(node)) {
    case FREESASA_NODE_STRUCTURE:
        *xml_node = freesasa_xml_structure(node);
        break;
    case FREESASA_NODE_CHAIN:
        *xml_node = freesasa_xml_chain(node);
        break;
    case FREESASA_NODE_RESIDUE:
        *xml_node = freesasa_xml_residue(node);
        break;
    case FREESASA_NODE_ATOM:
        *xml_node = freesasa_xml_atom(node);
        break;
    default:
        assert(0 && "Tree illegal");
    }
    if (!xml_node) return FREESASA_FAIL;

    // simplify?
    while (child) {
        if (freesasa_node2xml(&xml_child, child, exclude_type) == FREESASA_SUCCESS &&
            xml_child != NULL && 
            xmlAddChild(*xml_node, xml_child) == NULL) {
            fail_msg("");
            return FREESASA_FAIL;
        }
        child = freesasa_structure_node_next(child);
    }

    return FREESASA_SUCCESS;
}

static xmlNodePtr
parameters2xml(const freesasa_parameters *p)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "parameters");
    char buf[20];
    extern const char *freesasa_alg_names[];

    if (xml_node == NULL) {
        fail_msg("");
        return NULL;
    }

#ifdef HAVE_CONFIG_H
    if (xmlNewProp(xml_node, BAD_CAST "source", BAD_CAST PACKAGE_STRING) == NULL) {
        fail_msg("");
        goto cleanup;
    }
#endif

    if (xmlNewProp(xml_node, BAD_CAST "algorithm", freesasa_alg_names[p->alg]) == NULL) {
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

int
freesasa_write_xml(FILE *output,
                   const freesasa_structure_node *root,
                   const freesasa_parameters *parameters,
                   int options)
{
    freesasa_node_type exclude_type = FREESASA_NODE_NONE;
    xmlDocPtr doc = NULL; 
    xmlNodePtr xml_root = NULL, xml_structure = NULL, xml_param = NULL;
    xmlNsPtr ns = NULL;
    xmlBufferPtr buf = NULL;
    xmlTextWriterPtr writer = NULL;
    int ret = FREESASA_FAIL;

    if (parameters == NULL) parameters = &freesasa_default_parameters;
    if (options & FREESASA_OUTPUT_STRUCTURE) exclude_type = FREESASA_NODE_CHAIN;
    if (options & FREESASA_OUTPUT_CHAIN) exclude_type = FREESASA_NODE_RESIDUE;
    if (options & FREESASA_OUTPUT_RESIDUE) exclude_type = FREESASA_NODE_ATOM;

    doc = xmlNewDoc(BAD_CAST "1.0");
    if (doc == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_root = xmlNewNode(NULL, BAD_CAST "FreeSASAResult");
    if (xml_root == NULL) {
        fail_msg("");
        goto cleanup;
    }

    ns = xmlNewNs(xml_root, BAD_CAST "http://freesasa.github.io/", BAD_CAST "");
    buf = xmlBufferCreate();
    if (ns == NULL || buf == NULL) {
        fail_msg("");
        // this will later be stored with the doc, so can't be freed twice in cleanup
        xmlFreeNs(ns);
        goto cleanup;
    }

    xmlDocSetRootElement(doc, xml_root);

    writer = xmlNewTextWriterMemory(buf, 0);
    if (writer == NULL) {
        fail_msg("");
        goto cleanup;
    }

    xml_param = parameters2xml(parameters);
    if (xml_param == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_root, xml_param) == NULL) {
        fail_msg("");
        xmlFree(xml_param);
        goto cleanup;
    }

    if (xmlNewProp(xml_root, BAD_CAST "classifier",
                   BAD_CAST freesasa_structure_node_classified_by(root)) == NULL) {
        fail_msg("");
        goto cleanup;
    }
    
    if (xmlNewProp(xml_root, BAD_CAST "input",
                   BAD_CAST freesasa_structure_node_name(root)) == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlNewProp(xml_root, BAD_CAST "lengthUnit", BAD_CAST "Ångström") == NULL) {
        fail_msg("");
        goto cleanup;
    }

    if (freesasa_node2xml(&xml_structure, root, exclude_type) == FREESASA_FAIL) {
        fail_msg("");
        goto cleanup;
    }

    if (xmlAddChild(xml_root, xml_structure) == NULL) {
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
        fail_msg(buf->content);
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
    }
    ret = FREESASA_SUCCESS;
    
 cleanup:
    xmlFreeDoc(doc);
    xmlFreeTextWriter(writer);
    return ret;
}
