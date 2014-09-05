#include <stdlib.h>
#include <stdio.h>
#include "src/strmap.h"

int n_test;

void test(freesasa_strmap_t *rh) {
    printf("Test no: %d\n",++n_test);
    printf("ALA exists: %d\n",freesasa_strmap_exists(rh,"ALA"));
    printf("ARG exists: %d\n",freesasa_strmap_exists(rh,"ARG"));
    printf("ABCDEFG exists: %d\n",freesasa_strmap_exists(rh,"ABCDEFG"));
    printf("123 exists: %d\n",freesasa_strmap_exists(rh,"123"));
    printf("VAL exists: %d\n",freesasa_strmap_exists(rh,"VAL"));
    printf("size: %d\n",freesasa_strmap_size(rh));
    printf("\n");
}

int main(int argc, char **argv) {
    n_test = 0;
    freesasa_strmap_t *rh = freesasa_strmap_new();
    printf("Initialized\n");
    freesasa_strmap_real_add(rh,"ALA",1.0);
    freesasa_strmap_real_add(rh,"ALA",2.0);
    freesasa_strmap_set(rh,"ARG",NULL);
    freesasa_strmap_set(rh,"ABCDEFG",NULL);
    freesasa_strmap_set(rh,"123",NULL);
    freesasa_strmap_set(rh,"VAL",NULL);
    test(rh);
    freesasa_strmap_delete(rh,"ARG");
    freesasa_strmap_delete(rh,"123");
    test(rh);

    char **a = freesasa_strmap_keys(rh);
    for (int i = 0; i < freesasa_strmap_size(rh); ++i) {
        printf("%d %s\n",i,a[i]);
        free(a[i]);
    }
    free(a);
    printf("%f\n",freesasa_strmap_real_value(rh,"ALA"));

    freesasa_strmap_free(rh);
    return 0;
}
