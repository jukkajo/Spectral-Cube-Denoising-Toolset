#include "circular_right_shift.h"

#include <errno.h>
#include <stdlib.h>

#include "../Universal-Tools/universal_tools.h"

int circular_right_shift(double **input, size_t rows, size_t columns) {
    if (columns == 0U || validate_2d_array(input, rows) != 0) {
        errno = EINVAL;
        return -1;
    }
    if (columns == 1U) {
        return 0;
    }

    double *last_column = create_1d_array(rows);
    if (last_column == NULL) {
        return -1;
    }
    for (size_t row = 0U; row < rows; ++row) {
        last_column[row] = input[row][columns - 1U];
    }
    for (size_t column = columns - 1U; column > 0U; --column) {
        for (size_t row = 0U; row < rows; ++row) {
            input[row][column] = input[row][column - 1U];
        }
    }
    for (size_t row = 0U; row < rows; ++row) {
        input[row][0] = last_column[row];
    }
    free(last_column);
    return 0;
}
