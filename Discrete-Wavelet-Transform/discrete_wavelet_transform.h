#ifndef DISCRETE_WAVELET_TRANSFORM_H
#define DISCRETE_WAVELET_TRANSFORM_H

#include <stddef.h>

/*
 * One-level normalized Haar transform
 * -----------------------------------
 * For each adjacent pair (x[2i], x[2i + 1]):
 *
 *   approximation[i] = (x[2i] + x[2i + 1]) / sqrt(2)
 *   detail[i]        = (x[2i] - x[2i + 1]) / sqrt(2)
 *
 * An odd final sample is paired with a duplicate of itself. Consequently,
 * coefficient_length is ceil(original_length / 2), expressed without
 * overflowing as original_length / 2 + original_length % 2. The inverse
 * reconstructs every pair with
 *
 *   x[2i]     = (approximation[i] + detail[i]) / sqrt(2)
 *   x[2i + 1] = (approximation[i] - detail[i]) / sqrt(2)
 *
 * and omits the duplicated sample by returning original_length values.
 *
 * A successful forward transform returns a caller-owned result containing
 * two independent caller-owned coefficient arrays. The input is never
 * modified. Release the result and both arrays with haar_dwt_result_free().
 * The cleanup function accepts NULL.
 */
struct HaarDwtResult {
    size_t original_length;
    size_t coefficient_length;
    double *approximation;
    double *detail;
};

/*
 * Return NULL and set errno to EINVAL for a NULL input or zero input_length.
 * Allocation-size overflow uses EOVERFLOW; allocation failure uses ENOMEM.
 */
struct HaarDwtResult *haar_dwt_forward(
    const double *input,
    size_t input_length
);

/*
 * Validate all result metadata and coefficient pointers before reconstructing.
 * Invalid metadata returns NULL with errno set to EINVAL. Allocation-size
 * overflow uses EOVERFLOW; allocation failure uses ENOMEM. On success the
 * returned original_length-element array is caller-owned and must be freed
 * with free(). The result and its coefficients are not modified.
 */
double *haar_dwt_inverse(const struct HaarDwtResult *result);

void haar_dwt_result_free(struct HaarDwtResult *result);

/*
 * Reserved generic 3D forward-DWT entry point. No axis, subband layout, or
 * filter-bank phase convention is defined for this compatibility API. Valid
 * requests therefore return NULL with errno set to ENOTSUP; invalid requests
 * use EINVAL. No input memory is modified. Use haar_dwt_forward() for the
 * implemented one-dimensional Haar path.
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
