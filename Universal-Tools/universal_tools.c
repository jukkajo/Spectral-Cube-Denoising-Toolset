#include "universal_tools.h"

#include <errno.h>
#include <stdint.h>
#include <stdlib.h>

static int multiplication_overflows(size_t count, size_t element_size) {
    return element_size != 0U && count > SIZE_MAX / element_size;
}

static int invalid_dimensions(size_t rows, size_t columns) {
    if (rows == 0U || columns == 0U) {
        errno = EINVAL;
        return 1;
    }
    return 0;
}

int validate_2d_array(double *const *array, size_t rows) {
    if (array == NULL || rows == 0U) {
        errno = EINVAL;
        return -1;
    }

    for (size_t row = 0U; row < rows; ++row) {
        if (array[row] == NULL) {
            errno = EINVAL;
            return -1;
        }
    }
    return 0;
}

int validate_3d_array(double **const *array, size_t rows, size_t columns) {
    if (array == NULL || invalid_dimensions(rows, columns) != 0) {
        errno = EINVAL;
        return -1;
    }

    for (size_t row = 0U; row < rows; ++row) {
        if (array[row] == NULL) {
            errno = EINVAL;
            return -1;
        }
        for (size_t column = 0U; column < columns; ++column) {
            if (array[row][column] == NULL) {
                errno = EINVAL;
                return -1;
            }
        }
    }
    return 0;
}

void free_2d_array(double **array, size_t rows) {
    if (array == NULL) {
        return;
    }
    for (size_t row = 0U; row < rows; ++row) {
        free(array[row]);
    }
    free(array);
}

void free_3d_array(double ***array, size_t rows, size_t columns) {
    if (array == NULL) {
        return;
    }
    for (size_t row = 0U; row < rows; ++row) {
        if (array[row] != NULL) {
            for (size_t column = 0U; column < columns; ++column) {
                free(array[row][column]);
            }
        }
        free(array[row]);
    }
    free(array);
}

double *create_1d_array(size_t length) {
    if (length == 0U) {
        errno = EINVAL;
        return NULL;
    }
    if (multiplication_overflows(length, sizeof(double)) != 0) {
        errno = EOVERFLOW;
        return NULL;
    }

    double *array = calloc(length, sizeof(*array));
    if (array == NULL) {
        errno = ENOMEM;
    }
    return array;
}

double **create_2d_array(size_t rows, size_t columns) {
    if (invalid_dimensions(rows, columns) != 0) {
        return NULL;
    }
    if (multiplication_overflows(rows, sizeof(double *)) != 0 ||
        multiplication_overflows(columns, sizeof(double)) != 0) {
        errno = EOVERFLOW;
        return NULL;
    }

    double **array = calloc(rows, sizeof(*array));
    if (array == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    for (size_t row = 0U; row < rows; ++row) {
        array[row] = calloc(columns, sizeof(*array[row]));
        if (array[row] == NULL) {
            free_2d_array(array, rows);
            errno = ENOMEM;
            return NULL;
        }
    }
    return array;
}

double ***create_3d_zero_array(size_t rows, size_t columns, size_t pages) {
    if (invalid_dimensions(rows, columns) != 0 || pages == 0U) {
        errno = EINVAL;
        return NULL;
    }
    if (multiplication_overflows(rows, sizeof(double **)) != 0 ||
        multiplication_overflows(columns, sizeof(double *)) != 0 ||
        multiplication_overflows(pages, sizeof(double)) != 0) {
        errno = EOVERFLOW;
        return NULL;
    }

    double ***array = calloc(rows, sizeof(*array));
    if (array == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    for (size_t row = 0U; row < rows; ++row) {
        array[row] = calloc(columns, sizeof(*array[row]));
        if (array[row] == NULL) {
            free_3d_array(array, rows, columns);
            errno = ENOMEM;
            return NULL;
        }
        for (size_t column = 0U; column < columns; ++column) {
            array[row][column] = calloc(pages, sizeof(*array[row][column]));
            if (array[row][column] == NULL) {
                free_3d_array(array, rows, columns);
                errno = ENOMEM;
                return NULL;
            }
        }
    }
    return array;
}

double **duplicate_2d_array(double *const *source, size_t rows, size_t columns) {
    if (columns == 0U || validate_2d_array(source, rows) != 0) {
        errno = EINVAL;
        return NULL;
    }

    double **copy = create_2d_array(rows, columns);
    if (copy == NULL) {
        return NULL;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < columns; ++column) {
            copy[row][column] = source[row][column];
        }
    }
    return copy;
}

double **flatten_3d_to_2d(
    double **const *input,
    size_t rows,
    size_t columns,
    size_t pages
) {
    if (pages == 0U || validate_3d_array(input, rows, columns) != 0) {
        errno = EINVAL;
        return NULL;
    }
    if (columns > SIZE_MAX / rows) {
        errno = EOVERFLOW;
        return NULL;
    }

    const size_t flattened_rows = rows * columns;
    double **flattened = create_2d_array(flattened_rows, pages);
    if (flattened == NULL) {
        return NULL;
    }

    for (size_t column = 0U; column < columns; ++column) {
        for (size_t row = 0U; row < rows; ++row) {
            const size_t flat_row = column * rows + row;
            for (size_t page = 0U; page < pages; ++page) {
                flattened[flat_row][page] = input[row][column][page];
            }
        }
    }
    return flattened;
}

double ***reshape_2d_to_3d(
    double *const *input,
    size_t rows,
    size_t columns,
    size_t pages
) {
    if (invalid_dimensions(rows, columns) != 0 || pages == 0U ||
        columns > SIZE_MAX / rows) {
        errno = (columns > 0U && rows > 0U && columns > SIZE_MAX / rows)
            ? EOVERFLOW
            : EINVAL;
        return NULL;
    }

    const size_t flattened_rows = rows * columns;
    if (validate_2d_array(input, flattened_rows) != 0) {
        errno = EINVAL;
        return NULL;
    }

    double ***reshaped = create_3d_zero_array(rows, columns, pages);
    if (reshaped == NULL) {
        return NULL;
    }
    for (size_t column = 0U; column < columns; ++column) {
        for (size_t row = 0U; row < rows; ++row) {
            const size_t flat_row = column * rows + row;
            for (size_t page = 0U; page < pages; ++page) {
                reshaped[row][column][page] = input[flat_row][page];
            }
        }
    }
    return reshaped;
}

double **column_wise_2d_array_submatrix_duplication(
    double *const *source,
    size_t rows,
    size_t columns,
    size_t column_start,
    size_t column_end
) {
    if (columns == 0U || column_start >= column_end || column_end > columns ||
        validate_2d_array(source, rows) != 0) {
        errno = EINVAL;
        return NULL;
    }

    const size_t output_columns = column_end - column_start;
    double **submatrix = create_2d_array(rows, output_columns);
    if (submatrix == NULL) {
        return NULL;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < output_columns; ++column) {
            submatrix[row][column] = source[row][column_start + column];
        }
    }
    return submatrix;
}

double ***depth_wise_3d_submatrix_duplication(
    double **const *input,
    size_t rows,
    size_t columns,
    size_t pages,
    size_t page_start,
    size_t page_end
) {
    if (pages == 0U || page_start >= page_end || page_end > pages ||
        validate_3d_array(input, rows, columns) != 0) {
        errno = EINVAL;
        return NULL;
    }

    const size_t output_pages = page_end - page_start;
    double ***submatrix = create_3d_zero_array(rows, columns, output_pages);
    if (submatrix == NULL) {
        return NULL;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < columns; ++column) {
            for (size_t page = 0U; page < output_pages; ++page) {
                submatrix[row][column][page] =
                    input[row][column][page_start + page];
            }
        }
    }
    return submatrix;
}

int column_copying_worker(
    double *const *input,
    double **output,
    size_t rows,
    size_t columns,
    size_t output_columns,
    size_t step
) {
    if (columns == 0U || output_columns == 0U || step == 0U ||
        validate_2d_array(input, rows) != 0 ||
        validate_2d_array(output, rows) != 0) {
        errno = EINVAL;
        return -1;
    }
    if (columns - 1U > SIZE_MAX / step ||
        (columns - 1U) * step >= output_columns) {
        errno = EOVERFLOW;
        return -1;
    }

    for (size_t column = 0U; column < columns; ++column) {
        for (size_t row = 0U; row < rows; ++row) {
            output[row][column * step] = input[row][column];
        }
    }
    return 0;
}

double *reverse_1d_array(const double *array, size_t length) {
    if (array == NULL || length == 0U) {
        errno = EINVAL;
        return NULL;
    }

    double *reversed = create_1d_array(length);
    if (reversed == NULL) {
        return NULL;
    }
    for (size_t index = 0U; index < length; ++index) {
        reversed[index] = array[length - index - 1U];
    }
    return reversed;
}

double ***permute_3d(
    double **const *input,
    size_t rows,
    size_t columns,
    size_t pages,
    const size_t new_order[3]
) {
    if (pages == 0U || new_order == NULL ||
        validate_3d_array(input, rows, columns) != 0) {
        errno = EINVAL;
        return NULL;
    }

    unsigned int seen = 0U;
    for (size_t axis = 0U; axis < 3U; ++axis) {
        if (new_order[axis] < 1U || new_order[axis] > 3U) {
            errno = EINVAL;
            return NULL;
        }
        const unsigned int bit = 1U << (unsigned int)(new_order[axis] - 1U);
        if ((seen & bit) != 0U) {
            errno = EINVAL;
            return NULL;
        }
        seen |= bit;
    }

    const size_t input_dimensions[3] = {rows, columns, pages};
    const size_t output_rows = input_dimensions[new_order[0] - 1U];
    const size_t output_columns = input_dimensions[new_order[1] - 1U];
    const size_t output_pages = input_dimensions[new_order[2] - 1U];
    double ***output = create_3d_zero_array(
        output_rows,
        output_columns,
        output_pages
    );
    if (output == NULL) {
        return NULL;
    }

    for (size_t row = 0U; row < output_rows; ++row) {
        for (size_t column = 0U; column < output_columns; ++column) {
            for (size_t page = 0U; page < output_pages; ++page) {
                size_t source[3] = {0U, 0U, 0U};
                source[new_order[0] - 1U] = row;
                source[new_order[1] - 1U] = column;
                source[new_order[2] - 1U] = page;
                output[row][column][page] =
                    input[source[0]][source[1]][source[2]];
            }
        }
    }
    return output;
}
