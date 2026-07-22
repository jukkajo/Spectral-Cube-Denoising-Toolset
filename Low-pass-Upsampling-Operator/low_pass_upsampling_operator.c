#include "low_pass_upsampling_operator.h"

#include <errno.h>
#include <stdint.h>

#include "../Periodic-Convolution/periodic_convolution.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Upsampling-Operator/upsampling_operator.h"

double **low_pass_upsampling_operator(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t scale,
    size_t filter_length
) {
    if (columns == 0U || scale == 0U || filter == NULL ||
        filter_length == 0U || validate_2d_array(input, rows) != 0) {
        errno = EINVAL;
        return NULL;
    }
    if (columns > SIZE_MAX / scale) {
        errno = EOVERFLOW;
        return NULL;
    }

    double **upsampled = upsampling_operator(input, rows, columns, scale);
    if (upsampled == NULL) {
        return NULL;
    }
    const size_t output_columns = columns * scale;
    double **output = two_dim_convolution(
        filter,
        filter_length,
        upsampled,
        rows,
        output_columns
    );
    free_2d_array(upsampled, rows);
    return output;
}
