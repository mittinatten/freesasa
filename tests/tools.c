#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <dlfcn.h>
#include "tools.h"

static int n_fails = 0;
static int fail_after = 0;

void*
malloc(size_t s)
{
    if (fail_after > 0) {
        ++n_fails;
        if (n_fails >= fail_after) return NULL;
    }
    void *(*real_malloc)(size_t) = dlsym(RTLD_NEXT, "malloc");
    return real_malloc(s);
}


void*
realloc(void * ptr, size_t s)
{
    if (fail_after > 0) {
        ++n_fails;
        if (n_fails >= fail_after)
            return NULL;
    }
    void *(*real_realloc)(void *, size_t) = dlsym(RTLD_NEXT, "realloc");
    return real_realloc(ptr, s);
}

void*
strdup(const char *s)
{
    if (fail_after > 0) {
        ++n_fails;
        if (n_fails >= fail_after) {
            return NULL;
        }
    }
    void *(*real_strdup)(const char*) = dlsym(RTLD_NEXT, "strdup");
    return real_strdup(s);
}

void
set_fail_after(int limit) {
    if (limit < 0) limit = 0;
    fail_after = limit;
    n_fails = 0;
}

int float_eq(double a, double b, double tolerance) 
{
    if (fabs(a-b) < tolerance) return 1;
    printf("floats not equal: a = %f, b = %f, diff = %f, tolerance = %f\n",
           a, b, fabs(a-b), tolerance);
    fflush(stdout);
    return 0;
}
