#if HAVE_CONFIG_H
  #include <config.h>
#endif
#include <libxml/tree.h>
#include <libxml/xmlwriter.h>
#include <assert.h>
#include "freesasa_internal.h"

xmlNodePtr
freesasa_xml_structure(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "structure");
    xmlNewProp(xml_node, BAD_CAST "input", BAD_CAST freesasa_structure_node_name(node));
}

xmlNodePtr
freesasa_node2xml(const freesasa_structure_node *node, int exclude_type)
{
    assert(node);
    const freesasa_structure_node *child = freesasa_structure_node_children(node);
    xmlNodePtr xml_node = NULL, xml_child = NULL;

    switch (freesasa_structure_node_type(node)) {
    case FREESASA_NODE_STRUCTURE:
        xml_node = freesasa_xml_structure(node);
        break;
    case FREESASA_NODE_CHAIN:
        return NULL;
        break;
    case FREESASA_NODE_RESIDUE:
        return NULL;
        break;
    case FREESASA_NODE_ATOM:
        return NULL;
        break;
    default:
        assert(0 && "Tree illegal");
    }
    while (child) {
        xml_child = freesasa_node2xml(child, exclude_type);
        if (xml_child) {
            xmlAddChild(xml_node, xml_child);
        }
        child = freesasa_structure_node_next(child);
    }

    return xml_node;
}

int
freesasa_write_xml(FILE *output,
                   const freesasa_structure_node *root,
                   const freesasa_parameters *parameters,
                   int options)
{
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr xml_root = xmlNewNode(NULL, BAD_CAST "root");
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);

    xmlDocSetRootElement(doc, xml_root);
    xmlAddChild(xml_root, freesasa_node2xml(root,0));
    xmlNodeDump(buf, doc, xml_root, 0, 1);

    fprintf(output, "%s", (const char*) buf->content);
    
    return FREESASA_SUCCESS;
}
