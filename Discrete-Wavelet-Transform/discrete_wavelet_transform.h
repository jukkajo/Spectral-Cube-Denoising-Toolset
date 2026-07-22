#ifndef DISCRETE_WAVELET_TRANSFORM_H
#define DISCRETE_WAVELET_TRANSFORM_H

#include <stddef.h>

/*
 * Reserved forward-DWT entry point. Milestone 1 intentionally does not choose
 * an axis, subband layout, or filter-bank phase convention. Valid requests
 * therefore return NULL with errno set to ENOTSUP; invalid requests use
 * EINVAL. No input memory is modified.
 */
double ***discrete_wavelet_transform(
    double **const *input,
    const double *filter,
    size_t filter_length,
    size_t rows,
    size_t columns,
    size_t pages,
    size_t scale
);

#endif
