#ifndef RATIONAL_TRANSFER_FUNCTION_H
#define RATIONAL_TRANSFER_FUNCTION_H

#include <stddef.h>

/*
 * Causal row-wise FIR operation with zero extension:
 * output[r][n] = sum(filter[k] * input[r][n-k]) for k <= n.
 * All dimensions and filter_length must be positive. Returns a caller-owned
 * [rows][columns] array, or NULL with errno set on error. Input data and filter
 * coefficients are not modified.
 */
double **rational_transfer_function(
    double *const *input,
    size_t rows,
    size_t columns,
    const double *filter,
    size_t filter_length
);

#endif
