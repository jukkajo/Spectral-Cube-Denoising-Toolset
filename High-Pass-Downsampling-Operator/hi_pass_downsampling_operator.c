#include "hi_pass_downsampling_operator.h"

#include <errno.h>
#include <stdlib.h>

#include "../Circular-Left-Shift/circular_left_shift.h"
#include "../Mirror-Filter/mirror_filter.h"
#include "../Periodic-Convolution/periodic_convolution.h"
#include "../Universal-Tools/universal_tools.h"

double **hi_pass_downsampling_operator(
    double *const *input,
    const double *filter,
    size_t rows,
    size_t columns,
    size_t filter_length
) {
    if (columns == 0U || filter == NULL || filter_length == 0U ||
        validate_2d_array(input, rows) != 0) {
        errno = EINVAL;
        return NULL;
    }

    double *mirrored_filter = mirror_filter(filter, filter_length);
    if (mirrored_filter == NULL) {
        return NULL;
    }
    double **shifted_input = duplicate_2d_array(input, rows, columns);
    if (shifted_input == NULL) {
        free(mirrored_filter);
        return NULL;
    }
    if (circular_left_shift(shifted_input, rows, columns) != 0) {
        free_2d_array(shifted_input, rows);
        free(mirrored_filter);
        return NULL;
    }

    double **filtered = two_dim_convolution(
        mirrored_filter,
        filter_length,
        shifted_input,
        rows,
        columns
    );
    free_2d_array(shifted_input, rows);
    free(mirrored_filter);
    if (filtered == NULL) {
        return NULL;
    }

    const size_t output_columns = columns / 2U + columns % 2U;
    double **output = create_2d_array(rows, output_columns);
    if (output == NULL) {
        free_2d_array(filtered, rows);
        return NULL;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < output_columns; ++column) {
            output[row][column] = filtered[row][2U * column];
        }
    }
    free_2d_array(filtered, rows);
    return output;
}
