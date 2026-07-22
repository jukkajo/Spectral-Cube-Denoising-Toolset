#ifndef TIME_REVERSED_PERIODIC_CONVOLUTION_H
#define TIME_REVERSED_PERIODIC_CONVOLUTION_H

#include <stddef.h>

/*
 * Row-wise time-reversed periodic FIR operation:
 * output[r][n] = sum(filter[k] * input[r][(n+k) mod columns]).
 * All dimensions and filter_length must be positive. Filter lengths may be
 * shorter than, equal to, or longer than columns. Returns a caller-owned
 * [rows][columns] array, or NULL with errno set on error. Inputs are unchanged.
 */
double **time_reversed_periodic_convolution(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t filter_length
);

#endif
