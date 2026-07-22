#include "low_pass_downsampling_operator.h"

#include <errno.h>

#include "../Time-Reversed-Periodic-Convolution/time_reversed_periodic_convolution.h"
#include "../Universal-Tools/universal_tools.h"

double **low_pass_downsampling_operator(
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

    double **filtered = time_reversed_periodic_convolution(
        input,
        filter,
        rows,
        columns,
        filter_length
    );
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
