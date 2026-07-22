#ifndef HIGH_PASS_UPSAMPLING_OPERATOR_H
#define HIGH_PASS_UPSAMPLING_OPERATOR_H

#include <stddef.h>

/*
 * Upsample, circularly shift right, and periodically filter with a private
 * mirrored filter copy. Returns columns * scale columns without modifying
 * caller inputs. Dimensions, scale, and filter_length must be positive. NULL
 * is returned with errno set on error; successful output is caller-owned.
 */
double **high_pass_upsampling_operator(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t scale,
    size_t filter_length
);

#endif
