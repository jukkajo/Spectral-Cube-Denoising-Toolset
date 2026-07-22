#include "high_pass_upsampling_operator.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>

#include "../Circular-Right-Shift/circular_right_shift.h"
#include "../Mirror-Filter/mirror_filter.h"
#include "../Periodic-Convolution/periodic_convolution.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Upsampling-Operator/upsampling_operator.h"

double **high_pass_upsampling_operator(
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

    double *mirrored_filter = mirror_filter(filter, filter_length);
    if (mirrored_filter == NULL) {
        return NULL;
    }
    double **upsampled = upsampling_operator(input, rows, columns, scale);
    if (upsampled == NULL) {
        free(mirrored_filter);
        return NULL;
    }

    const size_t output_columns = columns * scale;
    if (circular_right_shift(upsampled, rows, output_columns) != 0) {
        free_2d_array(upsampled, rows);
        free(mirrored_filter);
        return NULL;
    }
    double **output = two_dim_convolution(
        mirrored_filter,
        filter_length,
        upsampled,
        rows,
        output_columns
    );
    free_2d_array(upsampled, rows);
    free(mirrored_filter);
    return output;
}
