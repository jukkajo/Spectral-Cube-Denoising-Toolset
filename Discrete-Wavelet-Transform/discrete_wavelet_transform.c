#include "discrete_wavelet_transform.h"

#include <errno.h>
#include <stdlib.h>

#include "../Universal-Tools/universal_tools.h"

#define HAAR_INV_SQRT_TWO 0.707106781186547524400844362104849039

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

    const size_t coefficient_length =
        input_length / 2U + input_length % 2U;
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
        result->original_length / 2U + result->original_length % 2U;
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
