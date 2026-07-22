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
 * Multilevel normalized Haar decomposition
 * ----------------------------------------
 * details[level], detail_lengths[level], and level_input_lengths[level]
 * describe the decomposition in forward order: index 0 is the first level
 * applied to the original signal, and index levels - 1 is the deepest level.
 * final_approximation is the approximation produced by the deepest level.
 * Every coefficient array is independently allocated.
 *
 * At each level, odd input lengths use the one-level duplicate-final-sample
 * rule documented above. level_input_lengths retains the unpadded length
 * needed to remove that level's reconstructed duplicate during inversion.
 */
struct HaarDwtDecomposition {
    size_t original_length;
    size_t levels;
    double *final_approximation;
    size_t final_approximation_length;
    double **details;
    size_t *detail_lengths;
    size_t *level_input_lengths;
};

/*
 * Return the largest meaningful multilevel count for input_length. A level
 * may start only with at least two samples; a one-sample approximation is
 * terminal and is not decomposed. The returned count is therefore zero for
 * input lengths zero and one. For larger inputs it is the number of repeated
 * ceil(length / 2) reductions needed to reach one sample.
 */
size_t haar_dwt_max_levels(size_t input_length);

/*
 * Decompose input through exactly levels levels by repeatedly invoking the
 * verified one-level Haar transform. input is not modified. Zero levels,
 * NULL input, zero input_length, and levels above haar_dwt_max_levels() return
 * NULL with errno set to EINVAL. Allocation arithmetic overflow uses
 * EOVERFLOW; allocation failure uses ENOMEM.
 *
 * On success the returned decomposition and all arrays it references are
 * caller-owned. Release the complete result with
 * haar_dwt_decomposition_free(), which accepts NULL and partially initialized
 * decomposition objects produced during allocation failure cleanup.
 */
struct HaarDwtDecomposition *haar_dwt_multilevel_forward(
    const double *input,
    size_t input_length,
    size_t levels
);

/*
 * Reconstruct the original signal by invoking the verified one-level inverse
 * from the deepest level to the first. All dimensions, per-level lengths, and
 * coefficient pointers are validated before reconstruction. Invalid metadata
 * returns NULL with errno set to EINVAL; allocation arithmetic overflow uses
 * EOVERFLOW; allocation failure uses ENOMEM. The decomposition is not
 * modified. The returned original_length-element array is caller-owned and
 * must be freed with free().
 */
double *haar_dwt_multilevel_inverse(
    const struct HaarDwtDecomposition *decomposition
);

void haar_dwt_decomposition_free(
    struct HaarDwtDecomposition *decomposition
);

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
