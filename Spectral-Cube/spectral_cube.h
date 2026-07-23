#ifndef SPECTRAL_CUBE_H
#define SPECTRAL_CUBE_H

#include <stddef.h>

#include "../Haar-Denoising/haar_denoising.h"

/*
 * Contiguous pixel-major spectral cube
 * ------------------------------------
 * data contains height * width * bands doubles in this order:
 *
 *   ((row * width) + column) * bands + band
 *
 * All band values belonging to one spatial pixel are therefore contiguous.
 * No internal row, page, or other pointer tree is exposed.
 */
typedef struct SpectralCube {
    size_t height;
    size_t width;
    size_t bands;
    double *data;
} SpectralCube;

/*
 * Store the checked height * width * bands count in element_count. Return 0
 * on success. Invalid cubes or a NULL output pointer return -1 with errno set
 * to EINVAL; dimension or byte-size arithmetic overflow uses EOVERFLOW.
 */
int spectral_cube_element_count(
    const SpectralCube *cube,
    size_t *element_count
);

/*
 * Create a cube with zero-initialized contiguous storage. All dimensions must
 * be nonzero. The returned cube is caller-owned and must be released with
 * spectral_cube_free(). Invalid dimensions use EINVAL, arithmetic overflow
 * uses EOVERFLOW, and allocation failure uses ENOMEM.
 */
SpectralCube *spectral_cube_create(
    size_t height,
    size_t width,
    size_t bands
);

/* Return an independent caller-owned copy of source. */
SpectralCube *spectral_cube_copy(const SpectralCube *source);

/* Release the cube and its contiguous data. Accepts NULL. */
void spectral_cube_free(SpectralCube *cube);

/*
 * Checked accessors using the documented pixel-major mapping. Return 0 on
 * success or -1 with errno set to EINVAL for NULL arguments, invalid cube
 * metadata, or an out-of-range row, column, or band. Metadata arithmetic
 * overflow uses EOVERFLOW.
 */
int spectral_cube_get(
    const SpectralCube *cube,
    size_t row,
    size_t column,
    size_t band,
    double *value
);

int spectral_cube_set(
    SpectralCube *cube,
    size_t row,
    size_t column,
    size_t band,
    double value
);

/*
 * Denoise each pixel's contiguous spectral vector independently with the
 * existing one-dimensional Haar denoising workflow. Only the spectral axis
 * is processed. Dimensions are preserved, source is never modified, and the
 * returned destination cube is independently caller-owned.
 *
 * bands must be at least two. levels, threshold, mode, and the optional
 * inclusive detail-level range follow haar_denoise_1d(). A NULL range selects
 * every detail level. Invalid arguments return NULL with errno set to EINVAL;
 * allocation arithmetic overflow uses EOVERFLOW and allocation failure uses
 * ENOMEM. All partial destination and temporary allocations are cleaned up on
 * failure.
 */
SpectralCube *spectral_cube_denoise_spectral(
    const SpectralCube *source,
    size_t levels,
    double threshold,
    enum HaarThresholdMode mode,
    const struct HaarDetailLevelRange *range
);

#endif
