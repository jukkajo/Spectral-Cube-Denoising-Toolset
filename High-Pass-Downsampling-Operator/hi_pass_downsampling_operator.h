#ifndef HIGH_PASS_DOWNSAMPLING_OPERATOR_H
#define HIGH_PASS_DOWNSAMPLING_OPERATOR_H

#include <stddef.h>

/*
 * Mirror a private filter copy, circularly shift a private input copy left,
 * apply periodic filtering, and retain columns 0, 2, 4, ... . The output has
 * columns / 2 + columns % 2 columns. Dimensions and filter_length must be
 * positive. NULL is returned with errno set on error. Caller inputs are not
 * modified.
 */
double **hi_pass_downsampling_operator(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t filter_length
);

#endif
