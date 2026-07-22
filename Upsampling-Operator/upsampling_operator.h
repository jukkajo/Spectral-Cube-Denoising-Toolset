#ifndef UPSAMPLING_OPERATOR_H
#define UPSAMPLING_OPERATOR_H

#include <stddef.h>

/*
 * Insert scale-1 zero columns after each input sample. The output contains
 * columns * scale columns and is caller-owned. Dimensions and scale must be
 * positive. Returns NULL with errno set to EINVAL, EOVERFLOW, or ENOMEM.
 */
double **upsampling_operator(
    double *const *input,
    size_t rows,
    size_t columns,
    size_t scale
);

#endif
