#ifndef HAAR_DENOISING_H
#define HAAR_DENOISING_H

#include <stddef.h>

#include "../Discrete-Wavelet-Transform/discrete_wavelet_transform.h"

enum HaarThresholdMode {
    HAAR_THRESHOLD_HARD = 0,
    HAAR_THRESHOLD_SOFT = 1
};

/*
 * Zero-based inclusive detail-level range in the decomposition's existing
 * forward ordering: level 0 is the first level and levels - 1 is deepest.
 */
struct HaarDetailLevelRange {
    size_t first_level;
    size_t last_level;
};

/*
 * Threshold coefficients in place using a finite, nonnegative threshold.
 *
 * Hard threshold:
 *   abs(value) < threshold  -> 0
 *   otherwise               -> value
 *
 * Soft threshold:
 *   sign(value) * max(abs(value) - threshold, 0)
 *
 * Equality is intentional: abs(value) == threshold remains unchanged for
 * hard thresholding and becomes zero for soft thresholding. Threshold zero
 * leaves coefficient values numerically unchanged in both modes.
 *
 * Return 0 on success. A NULL array, zero length, negative/NaN/infinite
 * threshold, or invalid mode returns -1 with errno set to EINVAL. Validation
 * completes before any coefficient is modified.
 */
int haar_threshold_coefficients(
    double *coefficients,
    size_t length,
    double threshold,
    enum HaarThresholdMode mode
);

/*
 * Threshold only decomposition detail coefficients in place. The final
 * approximation is never modified. A NULL range selects every detail level;
 * otherwise the inclusive range must satisfy
 * first_level <= last_level < decomposition->levels.
 *
 * Return 0 on success. Invalid decomposition metadata, threshold, mode, or
 * range returns -1 using the DWT errno conventions. All validation completes
 * before any detail coefficient is modified.
 */
int haar_threshold_detail_levels(
    struct HaarDwtDecomposition *decomposition,
    double threshold,
    enum HaarThresholdMode mode,
    const struct HaarDetailLevelRange *range
);

/*
 * Perform multilevel forward Haar decomposition, threshold the selected
 * detail levels, and reconstruct the one-dimensional signal. input is never
 * modified. A NULL range selects all levels. The threshold and level rules
 * are the same as above, and levels must satisfy the multilevel DWT contract;
 * in particular, zero levels and one-sample signals have no valid denoising
 * decomposition.
 *
 * On success the returned input_length-element array is caller-owned and must
 * be freed with free(). Invalid arguments return NULL with errno set to
 * EINVAL; allocation arithmetic overflow uses EOVERFLOW and allocation
 * failure uses ENOMEM. No threshold is selected automatically.
 */
double *haar_denoise_1d(
    const double *input,
    size_t input_length,
    size_t levels,
    double threshold,
    enum HaarThresholdMode mode,
    const struct HaarDetailLevelRange *range
);

#endif
