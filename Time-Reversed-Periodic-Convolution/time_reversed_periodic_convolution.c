#include "time_reversed_periodic_convolution.h"

#include <errno.h>

#include "../Universal-Tools/universal_tools.h"

static size_t add_periodically(size_t index, size_t amount, size_t length) {
    const size_t offset = amount % length;
    return (offset < length - index) ? index + offset : offset - (length - index);
}

double **time_reversed_periodic_convolution(
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

    double **output = create_2d_array(rows, columns);
    if (output == NULL) {
        return NULL;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < columns; ++column) {
            for (size_t tap = 0U; tap < filter_length; ++tap) {
                const size_t source = add_periodically(column, tap, columns);
                output[row][column] += input[row][source] * filter[tap];
            }
        }
    }
    return output;
}
