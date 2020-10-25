#include <check.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#if HAVE_CONFIG_H
#include <config.h>
#endif

extern Suite *pdb_suite();
extern Suite *classifier_suite();
extern Suite *coord_suite();
extern Suite *structure_suite();
extern Suite *sasa_suite();
extern Suite *nb_suite();
extern Suite *selector_suite();
extern Suite *result_node_suite();

#ifdef USE_JSON
extern Suite *json_suite();
#endif
#ifdef USE_XML
extern Suite *xml_suite();
#endif

int main(int argc, char **argv)
{
    mkdir("./tmp/", S_IRWXU);

    // Suites added in order of complexity
    SRunner *sr = srunner_create(pdb_suite());
    srunner_add_suite(sr, classifier_suite());
    srunner_add_suite(sr, coord_suite());
    srunner_add_suite(sr, structure_suite());
    srunner_add_suite(sr, sasa_suite());
    srunner_add_suite(sr, nb_suite());
    srunner_add_suite(sr, selector_suite());
    srunner_add_suite(sr, result_node_suite());
#if USE_JSON
    srunner_add_suite(sr, json_suite());
#endif
#if USE_XML
    srunner_add_suite(sr, xml_suite());
#endif
    srunner_run_all(sr, CK_VERBOSE);

    return (srunner_ntests_failed(sr) == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
