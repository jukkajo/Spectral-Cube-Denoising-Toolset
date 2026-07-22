#ifndef PERIODIC_CONVOLUTION_H
#define PERIODIC_CONVOLUTION_H

#include <stddef.h>

/*
 * Row-wise periodic FIR operation:
 * output[r][n] = sum(filter[k] * input[r][(n-k) mod columns]).
 * All dimensions and filter_length must be positive. Filter lengths may be
 * shorter than, equal to, or longer than columns. Returns a caller-owned
 * [rows][columns] array, or NULL with errno set on error. Inputs are unchanged.
 */
double **two_dim_convolution(
    const double *filter,
    size_t filter_length,
    double *const *input,
    size_t rows,
    size_t columns
);

#endif
