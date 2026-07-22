#include "discrete_wavelet_transform.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>

#include "../Universal-Tools/universal_tools.h"

#define HAAR_INV_SQRT_TWO 0.707106781186547524400844362104849039

static size_t half_rounded_up(size_t length) {
    return length / 2U + length % 2U;
}

size_t haar_dwt_max_levels(size_t input_length) {
    size_t levels = 0U;
    size_t length = input_length;

    while (length > 1U) {
        length = half_rounded_up(length);
        ++levels;
    }

    return levels;
}

static int add_allocation_size(
    size_t count,
    size_t element_size,
    size_t *total_size
) {
    if (element_size != 0U && count > SIZE_MAX / element_size) {
        errno = EOVERFLOW;
        return -1;
    }

    const size_t allocation_size = count * element_size;
    if (*total_size > SIZE_MAX - allocation_size) {
        errno = EOVERFLOW;
        return -1;
    }

    *total_size += allocation_size;
    return 0;
}

static int validate_multilevel_allocation_sizes(
    size_t input_length,
    size_t levels
) {
    size_t total_size = 0U;
    if (add_allocation_size(
            1U,
            sizeof(struct HaarDwtDecomposition),
            &total_size
        ) != 0 ||
        add_allocation_size(
            levels,
            sizeof(double *),
            &total_size
        ) != 0 ||
        add_allocation_size(
            levels,
            sizeof(size_t),
            &total_size
        ) != 0 ||
        add_allocation_size(
            levels,
            sizeof(size_t),
            &total_size
        ) != 0) {
        return -1;
    }

    size_t level_length = input_length;
    for (size_t level = 0U; level < levels; ++level) {
        level_length = half_rounded_up(level_length);
        if (add_allocation_size(
                level_length,
                sizeof(double),
                &total_size
            ) != 0) {
            return -1;
        }
    }

    return add_allocation_size(
        level_length,
        sizeof(double),
        &total_size
    );
}

void haar_dwt_result_free(struct HaarDwtResult *result) {
    if (result == NULL) {
        return;
    }

    free(result->approximation);
    free(result->detail);
    free(result);
}

struct HaarDwtResult *haar_dwt_forward(
    const double *input,
    size_t input_length
) {
    if (input == NULL || input_length == 0U) {
        errno = EINVAL;
        return NULL;
    }

    const size_t coefficient_length = half_rounded_up(input_length);
    struct HaarDwtResult *result = malloc(sizeof(*result));
    if (result == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    result->original_length = input_length;
    result->coefficient_length = coefficient_length;
    result->approximation = NULL;
    result->detail = NULL;

    result->approximation = create_1d_array(coefficient_length);
    if (result->approximation == NULL) {
        const int saved_errno = errno;
        haar_dwt_result_free(result);
        errno = saved_errno;
        return NULL;
    }

    result->detail = create_1d_array(coefficient_length);
    if (result->detail == NULL) {
        const int saved_errno = errno;
        haar_dwt_result_free(result);
        errno = saved_errno;
        return NULL;
    }

    for (size_t index = 0U; index < coefficient_length; ++index) {
        const size_t first_index = 2U * index;
        const size_t second_index = first_index + 1U;
        const double first = input[first_index];
        const double second = second_index < input_length
            ? input[second_index]
            : first;

        result->approximation[index] =
            (first + second) * HAAR_INV_SQRT_TWO;
        result->detail[index] =
            (first - second) * HAAR_INV_SQRT_TWO;
    }

    return result;
}

double *haar_dwt_inverse(const struct HaarDwtResult *result) {
    if (result == NULL || result->original_length == 0U ||
        result->coefficient_length == 0U ||
        result->approximation == NULL || result->detail == NULL) {
        errno = EINVAL;
        return NULL;
    }

    const size_t expected_coefficient_length =
        half_rounded_up(result->original_length);
    if (result->coefficient_length != expected_coefficient_length) {
        errno = EINVAL;
        return NULL;
    }

    double *output = create_1d_array(result->original_length);
    if (output == NULL) {
        return NULL;
    }

    for (size_t index = 0U; index < result->coefficient_length; ++index) {
        const size_t first_index = 2U * index;
        const size_t second_index = first_index + 1U;
        const double approximation = result->approximation[index];
        const double detail = result->detail[index];

        output[first_index] =
            (approximation + detail) * HAAR_INV_SQRT_TWO;
        if (second_index < result->original_length) {
            output[second_index] =
                (approximation - detail) * HAAR_INV_SQRT_TWO;
        }
    }

    return output;
}

void haar_dwt_decomposition_free(
    struct HaarDwtDecomposition *decomposition
) {
    if (decomposition == NULL) {
        return;
    }

    if (decomposition->details != NULL) {
        for (size_t level = 0U; level < decomposition->levels; ++level) {
            free(decomposition->details[level]);
        }
    }
    free(decomposition->details);
    free(decomposition->detail_lengths);
    free(decomposition->level_input_lengths);
    free(decomposition->final_approximation);
    free(decomposition);
}

static struct HaarDwtDecomposition *allocate_decomposition(
    size_t input_length,
    size_t levels
) {
    struct HaarDwtDecomposition *decomposition =
        malloc(sizeof(*decomposition));
    if (decomposition == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    decomposition->original_length = input_length;
    decomposition->levels = levels;
    decomposition->final_approximation = NULL;
    decomposition->final_approximation_length = 0U;
    decomposition->details = NULL;
    decomposition->detail_lengths = NULL;
    decomposition->level_input_lengths = NULL;

    decomposition->details = calloc(levels, sizeof(*decomposition->details));
    if (decomposition->details != NULL) {
        decomposition->detail_lengths =
            calloc(levels, sizeof(*decomposition->detail_lengths));
    }
    if (decomposition->detail_lengths != NULL) {
        decomposition->level_input_lengths =
            calloc(levels, sizeof(*decomposition->level_input_lengths));
    }
    if (decomposition->details == NULL ||
        decomposition->detail_lengths == NULL ||
        decomposition->level_input_lengths == NULL) {
        haar_dwt_decomposition_free(decomposition);
        errno = ENOMEM;
        return NULL;
    }

    return decomposition;
}

struct HaarDwtDecomposition *haar_dwt_multilevel_forward(
    const double *input,
    size_t input_length,
    size_t levels
) {
    if (input == NULL || input_length == 0U || levels == 0U ||
        levels > haar_dwt_max_levels(input_length)) {
        errno = EINVAL;
        return NULL;
    }
    if (validate_multilevel_allocation_sizes(input_length, levels) != 0) {
        return NULL;
    }

    struct HaarDwtDecomposition *decomposition =
        allocate_decomposition(input_length, levels);
    if (decomposition == NULL) {
        return NULL;
    }

    size_t level_input_length = input_length;
    for (size_t level = 0U; level < levels; ++level) {
        const double *level_input = level == 0U
            ? input
            : decomposition->final_approximation;
        struct HaarDwtResult *level_result =
            haar_dwt_forward(level_input, level_input_length);
        if (level_result == NULL) {
            const int saved_errno = errno;
            haar_dwt_decomposition_free(decomposition);
            errno = saved_errno;
            return NULL;
        }

        decomposition->level_input_lengths[level] = level_input_length;
        decomposition->detail_lengths[level] =
            level_result->coefficient_length;
        decomposition->details[level] = level_result->detail;
        level_result->detail = NULL;

        free(decomposition->final_approximation);
        decomposition->final_approximation = level_result->approximation;
        decomposition->final_approximation_length =
            level_result->coefficient_length;
        level_result->approximation = NULL;
        level_input_length = level_result->coefficient_length;
        haar_dwt_result_free(level_result);
    }

    return decomposition;
}

static int validate_decomposition(
    const struct HaarDwtDecomposition *decomposition
) {
    if (decomposition == NULL || decomposition->original_length < 2U ||
        decomposition->levels == 0U ||
        decomposition->levels >
            haar_dwt_max_levels(decomposition->original_length) ||
        decomposition->final_approximation == NULL ||
        decomposition->final_approximation_length == 0U ||
        decomposition->details == NULL ||
        decomposition->detail_lengths == NULL ||
        decomposition->level_input_lengths == NULL) {
        errno = EINVAL;
        return -1;
    }
    if (validate_multilevel_allocation_sizes(
            decomposition->original_length,
            decomposition->levels
        ) != 0) {
        return -1;
    }

    size_t expected_input_length = decomposition->original_length;
    for (size_t level = 0U; level < decomposition->levels; ++level) {
        const size_t expected_coefficient_length =
            half_rounded_up(expected_input_length);
        if (decomposition->level_input_lengths[level] !=
                expected_input_length ||
            decomposition->detail_lengths[level] !=
                expected_coefficient_length ||
            decomposition->details[level] == NULL) {
            errno = EINVAL;
            return -1;
        }
        expected_input_length = expected_coefficient_length;
    }

    if (decomposition->final_approximation_length != expected_input_length) {
        errno = EINVAL;
        return -1;
    }

    return 0;
}

double *haar_dwt_multilevel_inverse(
    const struct HaarDwtDecomposition *decomposition
) {
    if (validate_decomposition(decomposition) != 0) {
        return NULL;
    }

    double *current_approximation = decomposition->final_approximation;
    double *owned_approximation = NULL;
    for (size_t reverse_level = decomposition->levels;
         reverse_level > 0U;
         --reverse_level) {
        const size_t level = reverse_level - 1U;
        struct HaarDwtResult level_result = {
            decomposition->level_input_lengths[level],
            decomposition->detail_lengths[level],
            current_approximation,
            decomposition->details[level]
        };
        double *reconstructed = haar_dwt_inverse(&level_result);
        if (reconstructed == NULL) {
            const int saved_errno = errno;
            free(owned_approximation);
            errno = saved_errno;
            return NULL;
        }

        free(owned_approximation);
        owned_approximation = reconstructed;
        current_approximation = reconstructed;
    }

    return owned_approximation;
}

double ***discrete_wavelet_transform(
    double **const *input,
    const double *filter,
    size_t filter_length,
    size_t rows,
    size_t columns,
    size_t pages,
    size_t scale
) {
    if (pages == 0U || scale == 0U || filter == NULL || filter_length == 0U ||
        validate_3d_array(input, rows, columns) != 0) {
        errno = EINVAL;
        return NULL;
    }

    errno = ENOTSUP;
    return NULL;
}
