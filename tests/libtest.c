#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <check.h>

extern Suite* pdb_suite();
extern Suite* classifier_suite();
extern Suite* coord_suite();
extern Suite* structure_suite();
extern Suite* sasa_suite();
extern Suite* nb_suite();
extern Suite* selector_suite();

int main(int argc, char **argv) {
    mkdir("./tmp/",S_IRWXU);

    // Suites added in order of complexity
    SRunner *sr = srunner_create(pdb_suite());
    srunner_add_suite(sr,classifier_suite());
    srunner_add_suite(sr,coord_suite());
    srunner_add_suite(sr,structure_suite());
    srunner_add_suite(sr,sasa_suite());
    srunner_add_suite(sr,nb_suite());
    srunner_add_suite(sr,selector_suite());

    srunner_run_all(sr,CK_VERBOSE);

    return (srunner_ntests_failed(sr) == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
