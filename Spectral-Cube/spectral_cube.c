#include "spectral_cube.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int calculate_cube_sizes(
    size_t height,
    size_t width,
    size_t bands,
    size_t *element_count,
    size_t *byte_count
) {
    if (height == 0U || width == 0U || bands == 0U) {
        errno = EINVAL;
        return -1;
    }
    if (height > SIZE_MAX / width) {
        errno = EOVERFLOW;
        return -1;
    }

    const size_t pixel_count = height * width;
    if (pixel_count > SIZE_MAX / bands) {
        errno = EOVERFLOW;
        return -1;
    }

    const size_t total_elements = pixel_count * bands;
    if (total_elements > SIZE_MAX / sizeof(double)) {
        errno = EOVERFLOW;
        return -1;
    }

    *element_count = total_elements;
    *byte_count = total_elements * sizeof(double);
    return 0;
}

static int validate_cube(
    const SpectralCube *cube,
    size_t *element_count,
    size_t *byte_count
) {
    if (cube == NULL || cube->data == NULL) {
        errno = EINVAL;
        return -1;
    }
    return calculate_cube_sizes(
        cube->height,
        cube->width,
        cube->bands,
        element_count,
        byte_count
    );
}

int spectral_cube_element_count(
    const SpectralCube *cube,
    size_t *element_count
) {
    if (element_count == NULL) {
        errno = EINVAL;
        return -1;
    }

    size_t byte_count = 0U;
    return validate_cube(cube, element_count, &byte_count);
}

SpectralCube *spectral_cube_create(
    size_t height,
    size_t width,
    size_t bands
) {
    size_t element_count = 0U;
    size_t byte_count = 0U;
    if (calculate_cube_sizes(
            height,
            width,
            bands,
            &element_count,
            &byte_count
        ) != 0) {
        return NULL;
    }

    SpectralCube *cube = malloc(sizeof(*cube));
    if (cube == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    cube->height = height;
    cube->width = width;
    cube->bands = bands;
    cube->data = calloc(element_count, sizeof(*cube->data));
    if (cube->data == NULL) {
        free(cube);
        errno = ENOMEM;
        return NULL;
    }

    return cube;
}

SpectralCube *spectral_cube_copy(const SpectralCube *source) {
    size_t element_count = 0U;
    size_t byte_count = 0U;
    if (validate_cube(source, &element_count, &byte_count) != 0) {
        return NULL;
    }

    SpectralCube *copy = spectral_cube_create(
        source->height,
        source->width,
        source->bands
    );
    if (copy == NULL) {
        return NULL;
    }

    memcpy(copy->data, source->data, byte_count);
    return copy;
}

void spectral_cube_free(SpectralCube *cube) {
    if (cube == NULL) {
        return;
    }

    free(cube->data);
    free(cube);
}

static int checked_index(
    const SpectralCube *cube,
    size_t row,
    size_t column,
    size_t band,
    size_t *index
) {
    size_t element_count = 0U;
    size_t byte_count = 0U;
    if (validate_cube(cube, &element_count, &byte_count) != 0) {
        return -1;
    }
    if (row >= cube->height || column >= cube->width || band >= cube->bands) {
        errno = EINVAL;
        return -1;
    }

    *index = ((row * cube->width) + column) * cube->bands + band;
    return 0;
}

int spectral_cube_get(
    const SpectralCube *cube,
    size_t row,
    size_t column,
    size_t band,
    double *value
) {
    if (value == NULL) {
        errno = EINVAL;
        return -1;
    }

    size_t index = 0U;
    if (checked_index(cube, row, column, band, &index) != 0) {
        return -1;
    }

    *value = cube->data[index];
    return 0;
}

int spectral_cube_set(
    SpectralCube *cube,
    size_t row,
    size_t column,
    size_t band,
    double value
) {
    size_t index = 0U;
    if (checked_index(cube, row, column, band, &index) != 0) {
        return -1;
    }

    cube->data[index] = value;
    return 0;
}

SpectralCube *spectral_cube_denoise_spectral(
    const SpectralCube *source,
    size_t levels,
    double threshold,
    enum HaarThresholdMode mode,
    const struct HaarDetailLevelRange *range
) {
    size_t element_count = 0U;
    size_t byte_count = 0U;
    if (validate_cube(source, &element_count, &byte_count) != 0) {
        return NULL;
    }
    if (source->bands < 2U || levels == 0U ||
        levels > haar_dwt_max_levels(source->bands)) {
        errno = EINVAL;
        return NULL;
    }

    double *denoised_spectrum = haar_denoise_1d(
        source->data,
        source->bands,
        levels,
        threshold,
        mode,
        range
    );
    if (denoised_spectrum == NULL) {
        return NULL;
    }

    SpectralCube *destination = spectral_cube_create(
        source->height,
        source->width,
        source->bands
    );
    if (destination == NULL) {
        const int saved_errno = errno;
        free(denoised_spectrum);
        errno = saved_errno;
        return NULL;
    }

    const size_t spectrum_bytes = source->bands * sizeof(double);
    memcpy(destination->data, denoised_spectrum, spectrum_bytes);
    free(denoised_spectrum);

    const size_t pixel_count = element_count / source->bands;
    for (size_t pixel = 1U; pixel < pixel_count; ++pixel) {
        const size_t offset = pixel * source->bands;
        denoised_spectrum = haar_denoise_1d(
            source->data + offset,
            source->bands,
            levels,
            threshold,
            mode,
            range
        );
        if (denoised_spectrum == NULL) {
            const int saved_errno = errno;
            spectral_cube_free(destination);
            errno = saved_errno;
            return NULL;
        }

        memcpy(destination->data + offset, denoised_spectrum, spectrum_bytes);
        free(denoised_spectrum);
    }

    return destination;
}
