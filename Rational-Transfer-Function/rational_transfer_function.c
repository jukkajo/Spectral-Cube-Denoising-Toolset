#include "rational_transfer_function.h"

#include <errno.h>

#include "../Universal-Tools/universal_tools.h"

double **rational_transfer_function(
    double *const *input,
    size_t rows,
    size_t columns,
    const double *filter,
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
            const size_t terms = (filter_length < column + 1U)
                ? filter_length
                : column + 1U;
            for (size_t tap = 0U; tap < terms; ++tap) {
                output[row][column] += input[row][column - tap] * filter[tap];
            }
        }
    }
    return output;
}
