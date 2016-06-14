#if HAVE_CONFIG_H
  #include <config.h>
#endif
#include <libxml/tree.h>
#include <libxml/xmlwriter.h>
#include <string.h>
#include <assert.h>
#include "freesasa_internal.h"

static xmlNodePtr
xml_subarea(const freesasa_subarea *area, const char *name)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, name);
    char buf[20];

    sprintf(buf, "%f", area->total);
    xmlNewProp(xml_node, BAD_CAST "total", BAD_CAST buf);

    sprintf(buf, "%f", area->polar);
    xmlNewProp(xml_node, BAD_CAST "polar", BAD_CAST buf);

    sprintf(buf, "%f", area->apolar);
    xmlNewProp(xml_node, BAD_CAST "apolar", BAD_CAST buf);

    sprintf(buf, "%f", area->main_chain);
    xmlNewProp(xml_node, BAD_CAST "mainChain", BAD_CAST buf);

    sprintf(buf, "%f", area->side_chain);
    xmlNewProp(xml_node, BAD_CAST "sideChain", BAD_CAST buf);

    return xml_node;
}

xmlNodePtr
freesasa_xml_atom(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "atom");
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

    xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST trim_name);
    sprintf(buf, "%f", sasa);
    xmlNewProp(xml_node, BAD_CAST "area", BAD_CAST buf);
    sprintf(buf, "%s", is_polar ? "yes" : "no");
    xmlNewProp(xml_node, BAD_CAST "isPolar", BAD_CAST buf);
    sprintf(buf, "%s", is_bb ? "yes" : "no");
    xmlNewProp(xml_node, BAD_CAST "isMainChain", BAD_CAST buf);
    sprintf(buf, "%f", radius);
    xmlNewProp(xml_node, BAD_CAST "radius", BAD_CAST buf);

    return xml_node;
}

xmlNodePtr
freesasa_xml_residue(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = xmlNewNode(NULL, "residue");
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

    xmlNewProp(xml_node, BAD_CAST "name", BAD_CAST name);
    xmlNewProp(xml_node, BAD_CAST "number", BAD_CAST trim_number);
    xmlAddChild(xml_node, xml_subarea(abs, "area"));

    if (reference != NULL) {
        freesasa_residue_rel_subarea(&rel, abs, reference);
        xmlAddChild(xml_node, xml_subarea(&rel, "relativeArea"));
    }

    return xml_node;
}

xmlNodePtr
freesasa_xml_chain(const freesasa_structure_node *node)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "chain");
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    const char *name = freesasa_structure_node_name(node);
    char buf[20];
    int first, last;
    freesasa_structure_chain_residues(structure, name[0], &first, &last);

    xmlNewProp(xml_node, BAD_CAST "label", BAD_CAST freesasa_structure_node_name(node));
    sprintf(buf, "%d", last - first + 1);
    xmlNewProp(xml_node, BAD_CAST "nResidues", BAD_CAST buf);
    xmlAddChild(xml_node, xml_subarea(freesasa_structure_node_area(node), "area"));
    return xml_node;
}

xmlNodePtr
freesasa_xml_structure(const freesasa_structure_node *node)
{
    assert(node);
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "structure");
    const freesasa_structure *structure = freesasa_structure_node_structure(node);
    char buf[20];
    sprintf(buf, "%d", freesasa_structure_n_chains(structure));
    xmlNewProp(xml_node, BAD_CAST "nChains", BAD_CAST buf);
    xmlAddChild(xml_node, xml_subarea(freesasa_structure_node_area(node), "area"));
    return xml_node;
}

xmlNodePtr
freesasa_node2xml(const freesasa_structure_node *node, int exclude_type)
{
    assert(node);
    const freesasa_structure_node *child = freesasa_structure_node_children(node);
    xmlNodePtr xml_node = NULL, xml_child = NULL;

    if (freesasa_structure_node_type(node) == exclude_type) return NULL;;

    switch (freesasa_structure_node_type(node)) {
    case FREESASA_NODE_STRUCTURE:
        xml_node = freesasa_xml_structure(node);
        break;
    case FREESASA_NODE_CHAIN:
        xml_node = freesasa_xml_chain(node);
        break;
    case FREESASA_NODE_RESIDUE:
        xml_node = freesasa_xml_residue(node);
        break;
    case FREESASA_NODE_ATOM:
        xml_node = freesasa_xml_atom(node);
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

static xmlNodePtr
parameters2xml(const freesasa_parameters *p)
{
    xmlNodePtr xml_node = xmlNewNode(NULL, BAD_CAST "parameters");
    char buf[20];
    extern const char *freesasa_alg_names[];
#ifdef HAVE_CONFIG_H
    xmlNewProp(xml_node, BAD_CAST "source", BAD_CAST PACKAGE_STRING);
#endif
    xmlNewProp(xml_node, BAD_CAST "algorithm", freesasa_alg_names[p->alg]);
    sprintf(buf, "%f", p->probe_radius);
    xmlNewProp(xml_node, BAD_CAST "probeRadius", BAD_CAST buf);
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
    xmlNewProp(xml_node, BAD_CAST "resolution", BAD_CAST buf);
    return xml_node;
}

int
freesasa_write_xml(FILE *output,
                   const freesasa_structure_node *root,
                   const freesasa_parameters *parameters,
                   int options)
{
    freesasa_node_type exclude_type = FREESASA_NODE_NONE;
    xmlDocPtr doc = xmlNewDoc(BAD_CAST "1.0");
    xmlNodePtr xml_root = xmlNewNode(NULL, BAD_CAST "FreeSASAResult");
    xmlBufferPtr buf = xmlBufferCreate();
    xmlTextWriterPtr writer = xmlNewTextWriterMemory(buf, 0);
    if (parameters == NULL) parameters = &freesasa_default_parameters;

    if (options & FREESASA_OUTPUT_STRUCTURE) exclude_type = FREESASA_NODE_CHAIN;
    if (options & FREESASA_OUTPUT_CHAIN) exclude_type = FREESASA_NODE_RESIDUE;
    if (options & FREESASA_OUTPUT_RESIDUE) exclude_type = FREESASA_NODE_ATOM;

    xmlDocSetRootElement(doc, xml_root);
    xmlAddChild(xml_root, parameters2xml(parameters));
    xmlNewProp(xml_root, BAD_CAST "classifier",
               freesasa_structure_node_classified_by(root));
    xmlNewProp(xml_root, BAD_CAST "input",
               freesasa_structure_node_name(root));
    xmlNewProp(xml_root, BAD_CAST "lengthUnit", BAD_CAST "Ångström");
    xmlAddChild(xml_root, freesasa_node2xml(root, exclude_type));

    xmlTextWriterStartDocument(writer, "1.0", xmlGetCharEncodingName(XML_CHAR_ENCODING_UTF8), NULL);

    xmlTextWriterFlush(writer);
    xmlNodeDump(buf, doc, xml_root, 0, 1);
    xmlTextWriterEndDocument(writer);

    fprintf(output, "%s", (const char*) buf->content);

    xmlBufferFree(buf);
    
    return FREESASA_SUCCESS;
}
