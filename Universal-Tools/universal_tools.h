#ifndef UNIVERSAL_TOOLS_H
#define UNIVERSAL_TOOLS_H

#include <stddef.h>

/*
 * Array contract
 * --------------
 * Dimensions must be positive. Invalid arguments return NULL (allocating
 * functions) or -1 (status functions) and set errno. Every successful
 * allocating function returns caller-owned memory. The matching free helper
 * accepts NULL.
 *
 * 2D arrays are represented as rows independently allocated behind double **.
 * 3D arrays are represented as [rows][columns][pages] behind double ***.
 */

int validate_2d_array(double *const *array, size_t rows);
int validate_3d_array(double **const *array, size_t rows, size_t columns);

void free_2d_array(double **array, size_t rows);
void free_3d_array(double ***array, size_t rows, size_t columns);

double *create_1d_array(size_t length);
double **create_2d_array(size_t rows, size_t columns);
double ***create_3d_zero_array(size_t rows, size_t columns, size_t pages);
double **duplicate_2d_array(double *const *source, size_t rows, size_t columns);

/* Flatten [rows][columns][pages] to [rows * columns][pages]. */
double **flatten_3d_to_2d(
    double **const *input,
    size_t rows,
    size_t columns,
    size_t pages
);

/* Inverse of flatten_3d_to_2d for the supplied dimensions. */
double ***reshape_2d_to_3d(
    double *const *input,
    size_t rows,
    size_t columns,
    size_t pages
);

/* Copy the half-open column range [column_start, column_end). */
double **column_wise_2d_array_submatrix_duplication(
    double *const *source,
    size_t rows,
    size_t columns,
    size_t column_start,
    size_t column_end
);

/* Copy the half-open page range [page_start, page_end). */
double ***depth_wise_3d_submatrix_duplication(
    double **const *input,
    size_t rows,
    size_t columns,
    size_t pages,
    size_t page_start,
    size_t page_end
);

int column_copying_worker(
    double *const *input,
    double **output,
    size_t rows,
    size_t columns,
    size_t output_columns,
    size_t step
);

double *reverse_1d_array(const double *array, size_t length);

/*
 * new_order is a permutation of {1, 2, 3}; each value identifies the input
 * axis used for the corresponding output axis. The returned dimensions are
 * therefore dims[new_order[0]-1], dims[new_order[1]-1], dims[new_order[2]-1].
 */
double ***permute_3d(
    double **const *input,
    size_t rows,
    size_t columns,
    size_t pages,
    const size_t new_order[3]
);

#endif
