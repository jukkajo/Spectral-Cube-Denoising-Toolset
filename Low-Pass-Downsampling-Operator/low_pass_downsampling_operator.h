#ifndef LOW_PASS_DOWNSAMPLING_OPERATOR_H
#define LOW_PASS_DOWNSAMPLING_OPERATOR_H

#include <stddef.h>

/*
 * Apply time-reversed periodic filtering and retain columns 0, 2, 4, ... .
 * The caller-owned output has columns / 2 + columns % 2 columns. Dimensions
 * and filter_length must be positive. NULL is returned with errno set on error;
 * inputs are unchanged.
 */
double **low_pass_downsampling_operator(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t filter_length
);

#endif
