#include <stdio.h>
#include <math.h>
#include "tools.h"

#ifdef CHECK_MEMERR

static int n_fails = 0;
static int fail_after = 0;

extern void *__real_malloc(size_t size);
extern void *__real_realloc(void *ptr, size_t size);
extern void *__real_strdup(const char *s);

void*
__wrap_malloc(size_t s)
{
    if (fail_after > 0) {
        ++n_fails;
        if (n_fails >= fail_after) return NULL;
    }
    return __real_malloc(s);
}


void*
__wrap_realloc(void * ptr, size_t s)
{
    if (fail_after > 0) {
        ++n_fails;
        if (n_fails >= fail_after)
            return NULL;
    }
    return __real_realloc(ptr, s);
}

static void*
__wrap_strdup(const char *s)
{
    if (fail_after > 0) {
        ++n_fails;
        if (n_fails >= fail_after) return NULL;
    }
    return __real_strdup(s);
}

void
set_fail_after(int limit) {
    if (limit < 0) limit = 0;
    fail_after = limit;
    n_fails = 0;
}

#endif /* CHECK_MEMERR */

int float_eq(double a, double b, double tolerance) 
{
    if (fabs(a-b) < tolerance) return 1;
    printf("floats not equal: a = %f, b = %f, diff = %f, tolerance = %f\n",
           a, b, fabs(a-b), tolerance);
    fflush(stdout);
    return 0;
}
