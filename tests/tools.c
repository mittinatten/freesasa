#include <stdio.h>
#include <math.h>
#include "tools.h"


int float_eq(double a, double b, double tolerance) 
{
    if (fabs(a-b) < tolerance) return 1;
    printf("floats not equal: a = %f, b = %f, diff = %f, tolerance = %f\n",
           a, b, fabs(a-b), tolerance);
    fflush(stdout);
    return 0;
}
