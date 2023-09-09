#include <check.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#if HAVE_CONFIG_H
#include <config.h>
#endif

extern Suite *pdb_suite(void);
extern Suite *classifier_suite(void);
extern Suite *coord_suite(void);
extern Suite *structure_suite(void);
extern Suite *sasa_suite(void);
extern Suite *nb_suite(void);
extern Suite *selector_suite(void);
extern Suite *result_node_suite(void);

#ifdef USE_JSON
extern Suite *json_suite(void);
#endif
#ifdef USE_XML
extern Suite *xml_suite(void);
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
