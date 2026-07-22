#include "discrete_wavelet_transform.h"

#include <errno.h>

#include "../Universal-Tools/universal_tools.h"

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
