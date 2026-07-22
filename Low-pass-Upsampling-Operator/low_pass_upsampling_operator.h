#ifndef LOW_PASS_UPSAMPLING_OPERATOR_H
#define LOW_PASS_UPSAMPLING_OPERATOR_H

#include <stddef.h>

/*
 * Upsample, then periodically filter. Returns a caller-owned array containing
 * columns * scale columns. Dimensions, scale, and filter_length must be
 * positive. NULL is returned with errno set on error; inputs are unchanged.
 */
double **low_pass_upsampling_operator(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t scale,
    size_t filter_length
);

#endif
