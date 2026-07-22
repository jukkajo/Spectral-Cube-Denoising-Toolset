#include "upsampling_operator.h"

#include <errno.h>
#include <stdint.h>

#include "../Universal-Tools/universal_tools.h"

double **upsampling_operator(
    double *const *input,
    size_t rows,
    size_t columns,
    size_t scale
) {
    if (columns == 0U || scale == 0U || validate_2d_array(input, rows) != 0) {
        errno = EINVAL;
        return NULL;
    }
    if (columns > SIZE_MAX / scale) {
        errno = EOVERFLOW;
        return NULL;
    }

    const size_t output_columns = columns * scale;
    double **output = create_2d_array(rows, output_columns);
    if (output == NULL) {
        return NULL;
    }
    if (column_copying_worker(
            input,
            output,
            rows,
            columns,
            output_columns,
            scale
        ) != 0) {
        free_2d_array(output, rows);
        return NULL;
    }
    return output;
}
