#ifndef SRP_H
#define SRP_H

/** Prints the legal values for the number of points as a
    comma-separated list, ending with newline. */
void srp_print_n_opt(FILE*);

/** Returns an array of n test points (array has size 3*n). If n is
    not one of the legal values, an error message is printed and
    exit(1) is called. */
const double* srp_get_points(int n);

#endif
