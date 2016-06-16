#define _GNU_SOURCE
#include <freesasa_internal.h>
#include <check.h>
#include <libxml/tree.h>
#include <libxml/xmlwriter.h>
#include <dlfcn.h>
#include "tools.h"

#define fail_counter(err) if (fail_after > 0) { ++n_fails; if (n_fails >= fail_after) return (err); } 

static int n_fails;
static int fail_after;

static void
local_set_fail_after(int limit)
{
    if (limit < 0) limit = 0;
    fail_after = limit;
    n_fails = 0;
}

// mock functions
xmlDocPtr xmlNewDoc(const xmlChar *a)
{
    fail_counter(NULL);
    xmlDocPtr (*real_newdoc)(const xmlChar *) = dlsym(RTLD_NEXT, "xmlNewDoc");
    return real_newdoc(a);
}

xmlNodePtr xmlNewNode(xmlNsPtr a, const xmlChar *b)
{
    fail_counter(NULL);
    xmlNodePtr (*real_newnode)(xmlNsPtr, const xmlChar*) = dlsym(RTLD_NEXT, "xmlNewNode");
    return real_newnode(a, b);
}


xmlNsPtr xmlNewNs(xmlNodePtr a, const xmlChar *b, const xmlChar *c)
{
    fail_counter(NULL);
    xmlNsPtr (*real_newns)(xmlNodePtr, const xmlChar*, const xmlChar*) = dlsym(RTLD_NEXT, "xmlNewNs");
    return real_newns(a,b,c);
}

xmlBufferPtr xmlBufferCreate(void)
{
    fail_counter(NULL);
    xmlBufferPtr (*real_bc)() = dlsym(RTLD_NEXT, "xmlBufferCreate");
    return real_bc();
}

xmlTextWriterPtr xmlNewTextWriterMemory(xmlBufferPtr a, int b)
{
    fail_counter(NULL);
    xmlTextWriterPtr (*real_ntwm)(xmlBufferPtr, int) = dlsym(RTLD_NEXT, "xmlNewTextWriterMemory");
    return real_ntwm(a,b);
}

xmlNodePtr xmlAddChild(xmlNodePtr a, xmlNodePtr b) {
    fail_counter(NULL);
    xmlNodePtr (*real_add)(xmlNodePtr, xmlNodePtr) = dlsym(RTLD_NEXT, "xmlAddChild");
    return real_add(a, b);
}

xmlAttrPtr xmlNewProp(xmlNodePtr a, const xmlChar *b, const xmlChar *c)
{
    fail_counter(NULL);
    xmlAttrPtr (*real_np)(xmlNodePtr, const xmlChar*, const xmlChar*) = dlsym(RTLD_NEXT, "xmlNewProp");
    return real_np(a, b, c);
}

int
xmlTextWriterStartDocument(xmlTextWriterPtr a, const char *b, const char *c, const char *d)
{
    fail_counter(-1);
    int (*real_twsd)(xmlTextWriterPtr, const char*, const char*, const char*) = 
        dlsym(RTLD_NEXT, "xmlTextWriterStartDocument");
    return real_twsd(a, b, c, d);
}

int xmlTextWriterFlush(xmlTextWriterPtr a) {
    fail_counter(-1);
    int (*real_twf)(xmlTextWriterPtr) = dlsym(RTLD_NEXT, "xmlTextWriterFlush");
    return real_twf(a);
}

int xmlNodeDump(xmlBufferPtr a, xmlDocPtr b, xmlNodePtr c, int d, int e)
{
    fail_counter(0);
    int (*real_nd)(xmlBufferPtr, xmlDocPtr, xmlNodePtr, int, int) =
        dlsym(RTLD_NEXT, "xmlNodeDump");
    return real_nd(a, b, c, d, e);
}

int xmlTextWriterEndDocument(xmlTextWriterPtr a) {
    fail_counter(-1);
    int (*real_twed)(xmlTextWriterPtr) = dlsym(RTLD_NEXT, "xmlTextWriterEndDocument");
    return real_twed(a);
}

START_TEST (test_libxmlerr)
{
    FILE *pdb = fopen(DATADIR "1ubq.pdb", "r"), *devnull = fopen("/dev/null", "w");
    freesasa_structure *ubq =
        freesasa_structure_from_pdb(pdb, &freesasa_default_classifier, 0);
    fclose(pdb);
    freesasa_result *result = freesasa_calc_structure(ubq, NULL);
    freesasa_structure_node *root = freesasa_result2tree(result, ubq, NULL, "test");
    int ret;
    freesasa_set_verbosity(FREESASA_V_SILENT);
    for (int i = 1; i < 100; ++i) {
        local_set_fail_after(i);
        ret = freesasa_write_xml(devnull, root, NULL, FREESASA_OUTPUT_ATOM);
        local_set_fail_after(0);
        ck_assert_int_eq(ret, FREESASA_FAIL);
    }
    for (int i = 1; i < 29; ++i) {
        local_set_fail_after(i);
        ret = freesasa_write_xml(devnull, root, NULL, FREESASA_OUTPUT_STRUCTURE);
        local_set_fail_after(0);
        ck_assert_int_eq(ret, FREESASA_FAIL);
    }
    freesasa_set_verbosity(FREESASA_V_NORMAL);
}
END_TEST

Suite* xml_suite() {
    Suite *s = suite_create("XML");
    TCase *tc_core = tcase_create("Core");
    tcase_add_test(tc_core, test_libxmlerr);

    suite_add_tcase(s, tc_core);

    return s;
}
