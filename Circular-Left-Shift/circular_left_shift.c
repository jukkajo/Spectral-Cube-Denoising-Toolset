#include "circular_left_shift.h"

#include <errno.h>
#include <stdlib.h>

#include "../Universal-Tools/universal_tools.h"

int circular_left_shift(double **input, size_t rows, size_t columns) {
    if (columns == 0U || validate_2d_array(input, rows) != 0) {
        errno = EINVAL;
        return -1;
    }
    if (columns == 1U) {
        return 0;
    }

    double *first_column = create_1d_array(rows);
    if (first_column == NULL) {
        return -1;
    }
    for (size_t row = 0U; row < rows; ++row) {
        first_column[row] = input[row][0];
    }
    for (size_t column = 0U; column + 1U < columns; ++column) {
        for (size_t row = 0U; row < rows; ++row) {
            input[row][column] = input[row][column + 1U];
        }
    }
    for (size_t row = 0U; row < rows; ++row) {
        input[row][columns - 1U] = first_column[row];
    }
    free(first_column);
    return 0;
}
