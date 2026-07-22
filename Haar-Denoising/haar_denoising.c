#include "haar_denoising.h"

#include <errno.h>
#include <math.h>

static int validate_threshold(double threshold, enum HaarThresholdMode mode) {
    if (!isfinite(threshold) || threshold < 0.0 ||
        (mode != HAAR_THRESHOLD_HARD && mode != HAAR_THRESHOLD_SOFT)) {
        errno = EINVAL;
        return -1;
    }
    return 0;
}

static int validate_level_range(
    size_t levels,
    const struct HaarDetailLevelRange *range,
    size_t *first_level,
    size_t *last_level
) {
    if (levels == 0U) {
        errno = EINVAL;
        return -1;
    }

    if (range == NULL) {
        *first_level = 0U;
        *last_level = levels - 1U;
        return 0;
    }
    if (range->first_level > range->last_level ||
        range->last_level >= levels) {
        errno = EINVAL;
        return -1;
    }

    *first_level = range->first_level;
    *last_level = range->last_level;
    return 0;
}

static void apply_threshold(
    double *coefficients,
    size_t length,
    double threshold,
    enum HaarThresholdMode mode
) {
    for (size_t index = 0U; index < length; ++index) {
        const double value = coefficients[index];
        const double magnitude = fabs(value);

        if (mode == HAAR_THRESHOLD_HARD) {
            if (magnitude < threshold) {
                coefficients[index] = 0.0;
            }
        } else if (magnitude <= threshold) {
            coefficients[index] = 0.0;
        } else {
            coefficients[index] = copysign(magnitude - threshold, value);
        }
    }
}

int haar_threshold_coefficients(
    double *coefficients,
    size_t length,
    double threshold,
    enum HaarThresholdMode mode
) {
    if (coefficients == NULL || length == 0U ||
        validate_threshold(threshold, mode) != 0) {
        errno = EINVAL;
        return -1;
    }

    apply_threshold(coefficients, length, threshold, mode);
    return 0;
}

int haar_threshold_detail_levels(
    struct HaarDwtDecomposition *decomposition,
    double threshold,
    enum HaarThresholdMode mode,
    const struct HaarDetailLevelRange *range
) {
    size_t first_level = 0U;
    size_t last_level = 0U;
    if (validate_threshold(threshold, mode) != 0 ||
        haar_dwt_validate_decomposition(decomposition) != 0 ||
        validate_level_range(
            decomposition->levels,
            range,
            &first_level,
            &last_level
        ) != 0) {
        return -1;
    }

    for (size_t level = first_level; level <= last_level; ++level) {
        apply_threshold(
            decomposition->details[level],
            decomposition->detail_lengths[level],
            threshold,
            mode
        );
    }

    return 0;
}

double *haar_denoise_1d(
    const double *input,
    size_t input_length,
    size_t levels,
    double threshold,
    enum HaarThresholdMode mode,
    const struct HaarDetailLevelRange *range
) {
    size_t first_level = 0U;
    size_t last_level = 0U;
    if (input == NULL || input_length == 0U || levels == 0U ||
        levels > haar_dwt_max_levels(input_length) ||
        validate_threshold(threshold, mode) != 0 ||
        validate_level_range(
            levels,
            range,
            &first_level,
            &last_level
        ) != 0) {
        errno = EINVAL;
        return NULL;
    }

    struct HaarDwtDecomposition *decomposition =
        haar_dwt_multilevel_forward(input, input_length, levels);
    if (decomposition == NULL) {
        return NULL;
    }

    const struct HaarDetailLevelRange selected_range = {
        first_level,
        last_level
    };
    if (haar_threshold_detail_levels(
            decomposition,
            threshold,
            mode,
            &selected_range
        ) != 0) {
        const int saved_errno = errno;
        haar_dwt_decomposition_free(decomposition);
        errno = saved_errno;
        return NULL;
    }

    double *output = haar_dwt_multilevel_inverse(decomposition);
    if (output == NULL) {
        const int saved_errno = errno;
        haar_dwt_decomposition_free(decomposition);
        errno = saved_errno;
        return NULL;
    }

    haar_dwt_decomposition_free(decomposition);
    return output;
}
