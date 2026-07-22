#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "Circular-Left-Shift/circular_left_shift.h"
#include "Circular-Right-Shift/circular_right_shift.h"
#include "Discrete-Wavelet-Transform/discrete_wavelet_transform.h"
#include "High-Pass-Downsampling-Operator/hi_pass_downsampling_operator.h"
#include "High-Pass-Upsampling-Operator/high_pass_upsampling_operator.h"
#include "Low-Pass-Downsampling-Operator/low_pass_downsampling_operator.h"
#include "Low-pass-Upsampling-Operator/low_pass_upsampling_operator.h"
#include "Mirror-Filter/mirror_filter.h"
#include "Periodic-Convolution/periodic_convolution.h"
#include "Rational-Transfer-Function/rational_transfer_function.h"
#include "Time-Reversed-Periodic-Convolution/time_reversed_periodic_convolution.h"
#include "Universal-Tools/universal_tools.h"
#include "Upsampling-Operator/upsampling_operator.h"
#include "Wavelet-Coefficients/wavelet_coefficients.h"

static size_t checks_run = 0U;
static size_t failures = 0U;

static void check_condition(
    int condition,
    const char *expression,
    const char *file,
    int line
) {
    ++checks_run;
    if (condition == 0) {
        ++failures;
        fprintf(stderr, "%s:%d: check failed: %s\n", file, line, expression);
    }
}

#define CHECK(expression) \
    check_condition((expression), #expression, __FILE__, __LINE__)

static void check_close(double actual, double expected, const char *label) {
    ++checks_run;
    if (fabs(actual - expected) > 1.0e-12) {
        ++failures;
        fprintf(
            stderr,
            "%s: expected %.17g, got %.17g\n",
            label,
            expected,
            actual
        );
    }
}

static double **matrix_from_values(
    size_t rows,
    size_t columns,
    const double *values
) {
    double **matrix = create_2d_array(rows, columns);
    if (matrix == NULL) {
        return NULL;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < columns; ++column) {
            matrix[row][column] = values[row * columns + column];
        }
    }
    return matrix;
}

static void check_matrix_values(
    double *const *matrix,
    size_t rows,
    size_t columns,
    const double *expected,
    const char *label
) {
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < columns; ++column) {
            check_close(
                matrix[row][column],
                expected[row * columns + column],
                label
            );
        }
    }
}

static void test_allocation_and_copy_helpers(void) {
    double *vector = create_1d_array(3U);
    CHECK(vector != NULL);
    if (vector != NULL) {
        check_close(vector[0], 0.0, "zeroed 1D allocation");
        check_close(vector[1], 0.0, "zeroed 1D allocation");
        check_close(vector[2], 0.0, "zeroed 1D allocation");
    }
    free(vector);

    errno = 0;
    CHECK(create_1d_array(0U) == NULL);
    CHECK(errno == EINVAL);

    const double values[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double **matrix = matrix_from_values(2U, 3U, values);
    CHECK(matrix != NULL);
    if (matrix == NULL) {
        return;
    }
    CHECK(validate_2d_array(matrix, 2U) == 0);

    double **copy = duplicate_2d_array(matrix, 2U, 3U);
    CHECK(copy != NULL);
    if (copy != NULL) {
        check_matrix_values(copy, 2U, 3U, values, "2D duplicate");
        copy[0][0] = 99.0;
        check_close(matrix[0][0], 1.0, "2D duplicate ownership");
    }
    free_2d_array(copy, 2U);

    double **columns = column_wise_2d_array_submatrix_duplication(
        matrix,
        2U,
        3U,
        1U,
        3U
    );
    const double expected_columns[] = {2.0, 3.0, 5.0, 6.0};
    CHECK(columns != NULL);
    if (columns != NULL) {
        check_matrix_values(
            columns,
            2U,
            2U,
            expected_columns,
            "column submatrix"
        );
    }
    free_2d_array(columns, 2U);

    double **worker_output = create_2d_array(2U, 5U);
    CHECK(worker_output != NULL);
    if (worker_output != NULL) {
        CHECK(column_copying_worker(matrix, worker_output, 2U, 3U, 5U, 2U) == 0);
        const double expected_worker[] = {
            1.0, 0.0, 2.0, 0.0, 3.0,
            4.0, 0.0, 5.0, 0.0, 6.0
        };
        check_matrix_values(
            worker_output,
            2U,
            5U,
            expected_worker,
            "column worker"
        );
        errno = 0;
        CHECK(column_copying_worker(matrix, worker_output, 2U, 3U, 4U, 2U) == -1);
        CHECK(errno == EOVERFLOW);
    }
    free_2d_array(worker_output, 2U);

    double *reversed = reverse_1d_array(values, 6U);
    CHECK(reversed != NULL);
    if (reversed != NULL) {
        for (size_t index = 0U; index < 6U; ++index) {
            check_close(reversed[index], values[5U - index], "reverse 1D");
        }
    }
    free(reversed);

    double *broken_rows[2] = {matrix[0], NULL};
    errno = 0;
    CHECK(validate_2d_array(broken_rows, 2U) == -1);
    CHECK(errno == EINVAL);

    errno = 0;
    CHECK(create_2d_array(0U, 2U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(create_2d_array(1U, SIZE_MAX) == NULL);
    CHECK(errno == EOVERFLOW);
    errno = 0;
    CHECK(column_wise_2d_array_submatrix_duplication(
        matrix,
        2U,
        3U,
        2U,
        2U
    ) == NULL);
    CHECK(errno == EINVAL);

    free_2d_array(matrix, 2U);
    free_2d_array(NULL, 0U);
    free_3d_array(NULL, 0U, 0U);
}

static void test_3d_helpers(void) {
    const size_t rows = 2U;
    const size_t columns = 2U;
    const size_t pages = 3U;
    double ***cube = create_3d_zero_array(rows, columns, pages);
    CHECK(cube != NULL);
    if (cube == NULL) {
        return;
    }
    for (size_t row = 0U; row < rows; ++row) {
        for (size_t column = 0U; column < columns; ++column) {
            for (size_t page = 0U; page < pages; ++page) {
                cube[row][column][page] =
                    (double)(100U * row + 10U * column + page);
            }
        }
    }
    CHECK(validate_3d_array(cube, rows, columns) == 0);

    double **flat = flatten_3d_to_2d(cube, rows, columns, pages);
    CHECK(flat != NULL);
    if (flat != NULL) {
        for (size_t column = 0U; column < columns; ++column) {
            for (size_t row = 0U; row < rows; ++row) {
                const size_t flat_row = column * rows + row;
                for (size_t page = 0U; page < pages; ++page) {
                    check_close(
                        flat[flat_row][page],
                        cube[row][column][page],
                        "flatten 3D"
                    );
                }
            }
        }

        double ***reshaped = reshape_2d_to_3d(flat, rows, columns, pages);
        CHECK(reshaped != NULL);
        if (reshaped != NULL) {
            for (size_t row = 0U; row < rows; ++row) {
                for (size_t column = 0U; column < columns; ++column) {
                    for (size_t page = 0U; page < pages; ++page) {
                        check_close(
                            reshaped[row][column][page],
                            cube[row][column][page],
                            "reshape 2D to 3D"
                        );
                    }
                }
            }
        }
        free_3d_array(reshaped, rows, columns);
    }
    free_2d_array(flat, rows * columns);

    double ***depth = depth_wise_3d_submatrix_duplication(
        cube,
        rows,
        columns,
        pages,
        1U,
        3U
    );
    CHECK(depth != NULL);
    if (depth != NULL) {
        for (size_t row = 0U; row < rows; ++row) {
            for (size_t column = 0U; column < columns; ++column) {
                check_close(
                    depth[row][column][0],
                    cube[row][column][1],
                    "depth slice first"
                );
                check_close(
                    depth[row][column][1],
                    cube[row][column][2],
                    "depth slice second"
                );
            }
        }
    }
    free_3d_array(depth, rows, columns);

    const size_t orders[6][3] = {
        {1U, 2U, 3U}, {1U, 3U, 2U}, {2U, 1U, 3U},
        {2U, 3U, 1U}, {3U, 1U, 2U}, {3U, 2U, 1U}
    };
    const size_t dimensions[3] = {rows, columns, pages};
    for (size_t permutation = 0U; permutation < 6U; ++permutation) {
        const size_t *order = orders[permutation];
        const size_t output_rows = dimensions[order[0] - 1U];
        const size_t output_columns = dimensions[order[1] - 1U];
        const size_t output_pages = dimensions[order[2] - 1U];
        double ***permuted = permute_3d(
            cube,
            rows,
            columns,
            pages,
            order
        );
        CHECK(permuted != NULL);
        if (permuted != NULL) {
            for (size_t row = 0U; row < output_rows; ++row) {
                for (size_t column = 0U; column < output_columns; ++column) {
                    for (size_t page = 0U; page < output_pages; ++page) {
                        size_t source[3] = {0U, 0U, 0U};
                        source[order[0] - 1U] = row;
                        source[order[1] - 1U] = column;
                        source[order[2] - 1U] = page;
                        check_close(
                            permuted[row][column][page],
                            cube[source[0]][source[1]][source[2]],
                            "3D permutation"
                        );
                    }
                }
            }
        }
        free_3d_array(permuted, output_rows, output_columns);
    }

    const size_t invalid_order[3] = {1U, 1U, 3U};
    errno = 0;
    CHECK(permute_3d(cube, rows, columns, pages, invalid_order) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(depth_wise_3d_submatrix_duplication(
        cube,
        rows,
        columns,
        pages,
        2U,
        2U
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(create_3d_zero_array(SIZE_MAX, 1U, 1U) == NULL);
    CHECK(errno == EOVERFLOW);

    free_3d_array(cube, rows, columns);
}

static void test_circular_shifts(void) {
    const double original[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const double shifted_left[] = {2.0, 3.0, 1.0, 5.0, 6.0, 4.0};
    double **matrix = matrix_from_values(2U, 3U, original);
    CHECK(matrix != NULL);
    if (matrix == NULL) {
        return;
    }
    CHECK(circular_left_shift(matrix, 2U, 3U) == 0);
    check_matrix_values(matrix, 2U, 3U, shifted_left, "circular left shift");
    CHECK(circular_right_shift(matrix, 2U, 3U) == 0);
    check_matrix_values(matrix, 2U, 3U, original, "circular right shift");
    free_2d_array(matrix, 2U);

    const double singleton_value[] = {7.0};
    double **singleton = matrix_from_values(1U, 1U, singleton_value);
    CHECK(singleton != NULL);
    if (singleton != NULL) {
        CHECK(circular_left_shift(singleton, 1U, 1U) == 0);
        CHECK(circular_right_shift(singleton, 1U, 1U) == 0);
        check_close(singleton[0][0], 7.0, "singleton circular shift");
    }
    free_2d_array(singleton, 1U);

    errno = 0;
    CHECK(circular_left_shift(NULL, 0U, 0U) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(circular_right_shift(NULL, 0U, 0U) == -1);
    CHECK(errno == EINVAL);
}

static void test_filtering_and_convolution(void) {
    const double values[] = {1.0, 2.0, 3.0, 4.0};
    const double two_tap[] = {1.0, 1.0};
    double **input = matrix_from_values(1U, 4U, values);
    CHECK(input != NULL);
    if (input == NULL) {
        return;
    }

    double **causal = rational_transfer_function(input, 1U, 4U, two_tap, 2U);
    const double causal_expected[] = {1.0, 3.0, 5.0, 7.0};
    CHECK(causal != NULL);
    if (causal != NULL) {
        check_matrix_values(causal, 1U, 4U, causal_expected, "causal FIR");
    }
    free_2d_array(causal, 1U);

    double **periodic = two_dim_convolution(two_tap, 2U, input, 1U, 4U);
    const double periodic_expected[] = {5.0, 3.0, 5.0, 7.0};
    CHECK(periodic != NULL);
    if (periodic != NULL) {
        check_matrix_values(periodic, 1U, 4U, periodic_expected, "periodic FIR");
    }
    free_2d_array(periodic, 1U);

    double **time_reversed = time_reversed_periodic_convolution(
        input,
        two_tap,
        1U,
        4U,
        2U
    );
    const double reversed_expected[] = {3.0, 5.0, 7.0, 5.0};
    CHECK(time_reversed != NULL);
    if (time_reversed != NULL) {
        check_matrix_values(
            time_reversed,
            1U,
            4U,
            reversed_expected,
            "time-reversed periodic FIR"
        );
    }
    free_2d_array(time_reversed, 1U);

    const double equal_length_filter[] = {1.0, 0.0, 0.0, 0.0};
    periodic = two_dim_convolution(equal_length_filter, 4U, input, 1U, 4U);
    CHECK(periodic != NULL);
    if (periodic != NULL) {
        check_matrix_values(periodic, 1U, 4U, values, "equal-length filter");
    }
    free_2d_array(periodic, 1U);

    const double long_filter[] = {1.0, 0.0, 0.0, 0.0, 1.0};
    const double doubled[] = {2.0, 4.0, 6.0, 8.0};
    periodic = two_dim_convolution(long_filter, 5U, input, 1U, 4U);
    CHECK(periodic != NULL);
    if (periodic != NULL) {
        check_matrix_values(periodic, 1U, 4U, doubled, "long periodic filter");
    }
    free_2d_array(periodic, 1U);

    time_reversed = time_reversed_periodic_convolution(
        input,
        long_filter,
        1U,
        4U,
        5U
    );
    CHECK(time_reversed != NULL);
    if (time_reversed != NULL) {
        check_matrix_values(
            time_reversed,
            1U,
            4U,
            doubled,
            "long time-reversed filter"
        );
    }
    free_2d_array(time_reversed, 1U);

    const double singleton_value[] = {2.0};
    const double singleton_filter[] = {1.0, 2.0, 3.0};
    double **singleton = matrix_from_values(1U, 1U, singleton_value);
    CHECK(singleton != NULL);
    if (singleton != NULL) {
        periodic = two_dim_convolution(singleton_filter, 3U, singleton, 1U, 1U);
        CHECK(periodic != NULL);
        if (periodic != NULL) {
            check_close(periodic[0][0], 12.0, "singleton long periodic filter");
        }
        free_2d_array(periodic, 1U);
    }
    free_2d_array(singleton, 1U);

    check_matrix_values(input, 1U, 4U, values, "convolution input ownership");
    check_close(two_tap[0], 1.0, "convolution filter ownership");
    check_close(two_tap[1], 1.0, "convolution filter ownership");

    errno = 0;
    CHECK(two_dim_convolution(two_tap, 2U, NULL, 0U, 0U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(time_reversed_periodic_convolution(NULL, two_tap, 0U, 0U, 2U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(rational_transfer_function(input, 1U, 4U, NULL, 0U) == NULL);
    CHECK(errno == EINVAL);

    free_2d_array(input, 1U);
}

static void test_mirror_and_sampling(void) {
    const double filter[] = {1.0, 2.0, 3.0, 4.0};
    const double mirrored_expected[] = {1.0, -2.0, 3.0, -4.0};
    double *mirrored = mirror_filter(filter, 4U);
    CHECK(mirrored != NULL);
    if (mirrored != NULL) {
        for (size_t index = 0U; index < 4U; ++index) {
            check_close(mirrored[index], mirrored_expected[index], "mirror filter");
            check_close(filter[index], (double)(index + 1U), "mirror input ownership");
        }
        double *restored = mirror_filter(mirrored, 4U);
        CHECK(restored != NULL);
        if (restored != NULL) {
            for (size_t index = 0U; index < 4U; ++index) {
                check_close(restored[index], filter[index], "double mirror");
            }
        }
        free(restored);
    }
    free(mirrored);
    errno = 0;
    CHECK(mirror_filter(NULL, 0U) == NULL);
    CHECK(errno == EINVAL);

    const double odd_values[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    const double identity_filter[] = {1.0};
    double **odd = matrix_from_values(1U, 5U, odd_values);
    CHECK(odd != NULL);
    if (odd == NULL) {
        return;
    }

    double **upsampled = upsampling_operator(odd, 1U, 5U, 2U);
    const double upsampled_expected[] = {
        1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 5.0, 0.0
    };
    CHECK(upsampled != NULL);
    if (upsampled != NULL) {
        check_matrix_values(
            upsampled,
            1U,
            10U,
            upsampled_expected,
            "upsampling"
        );
    }
    free_2d_array(upsampled, 1U);

    double **low_down = low_pass_downsampling_operator(
        odd,
        identity_filter,
        1U,
        5U,
        1U
    );
    const double low_down_expected[] = {1.0, 3.0, 5.0};
    CHECK(low_down != NULL);
    if (low_down != NULL) {
        check_matrix_values(
            low_down,
            1U,
            3U,
            low_down_expected,
            "odd low-pass downsampling"
        );
    }
    free_2d_array(low_down, 1U);

    double **high_down = hi_pass_downsampling_operator(
        odd,
        identity_filter,
        1U,
        5U,
        1U
    );
    const double high_down_expected[] = {2.0, 4.0, 1.0};
    CHECK(high_down != NULL);
    if (high_down != NULL) {
        check_matrix_values(
            high_down,
            1U,
            3U,
            high_down_expected,
            "odd high-pass downsampling"
        );
    }
    free_2d_array(high_down, 1U);
    check_matrix_values(odd, 1U, 5U, odd_values, "downsampling input ownership");
    check_close(identity_filter[0], 1.0, "downsampling filter ownership");

    const double even_values[] = {1.0, 2.0, 3.0, 4.0};
    double **even = matrix_from_values(1U, 4U, even_values);
    CHECK(even != NULL);
    if (even != NULL) {
        low_down = low_pass_downsampling_operator(
            even,
            identity_filter,
            1U,
            4U,
            1U
        );
        const double even_expected[] = {1.0, 3.0};
        CHECK(low_down != NULL);
        if (low_down != NULL) {
            check_matrix_values(
                low_down,
                1U,
                2U,
                even_expected,
                "even low-pass downsampling"
            );
        }
        free_2d_array(low_down, 1U);
    }
    free_2d_array(even, 1U);

    const double pair_values[] = {1.0, 2.0};
    double **pair = matrix_from_values(1U, 2U, pair_values);
    CHECK(pair != NULL);
    if (pair != NULL) {
        double **low_up = low_pass_upsampling_operator(
            pair,
            identity_filter,
            1U,
            2U,
            2U,
            1U
        );
        const double low_up_expected[] = {1.0, 0.0, 2.0, 0.0};
        CHECK(low_up != NULL);
        if (low_up != NULL) {
            check_matrix_values(
                low_up,
                1U,
                4U,
                low_up_expected,
                "low-pass upsampling"
            );
        }
        free_2d_array(low_up, 1U);

        double **high_up = high_pass_upsampling_operator(
            pair,
            identity_filter,
            1U,
            2U,
            2U,
            1U
        );
        const double high_up_expected[] = {0.0, 1.0, 0.0, 2.0};
        CHECK(high_up != NULL);
        if (high_up != NULL) {
            check_matrix_values(
                high_up,
                1U,
                4U,
                high_up_expected,
                "high-pass upsampling"
            );
        }
        free_2d_array(high_up, 1U);
        check_matrix_values(pair, 1U, 2U, pair_values, "upsampling input ownership");
    }
    free_2d_array(pair, 1U);

    const double singleton_value[] = {7.0};
    double **singleton = matrix_from_values(1U, 1U, singleton_value);
    CHECK(singleton != NULL);
    if (singleton != NULL) {
        low_down = low_pass_downsampling_operator(
            singleton,
            identity_filter,
            1U,
            1U,
            1U
        );
        high_down = hi_pass_downsampling_operator(
            singleton,
            identity_filter,
            1U,
            1U,
            1U
        );
        CHECK(low_down != NULL);
        CHECK(high_down != NULL);
        if (low_down != NULL) {
            check_close(low_down[0][0], 7.0, "singleton low downsampling");
        }
        if (high_down != NULL) {
            check_close(high_down[0][0], 7.0, "singleton high downsampling");
        }
        free_2d_array(low_down, 1U);
        free_2d_array(high_down, 1U);
    }
    free_2d_array(singleton, 1U);

    errno = 0;
    CHECK(upsampling_operator(odd, 1U, 5U, 0U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(upsampling_operator(odd, 1U, SIZE_MAX, 2U) == NULL);
    CHECK(errno == EOVERFLOW);
    errno = 0;
    CHECK(low_pass_downsampling_operator(NULL, identity_filter, 0U, 0U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(hi_pass_downsampling_operator(NULL, identity_filter, 0U, 0U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(low_pass_upsampling_operator(NULL, identity_filter, 0U, 0U, 2U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(high_pass_upsampling_operator(NULL, identity_filter, 0U, 0U, 2U, 1U) == NULL);
    CHECK(errno == EINVAL);

    free_2d_array(odd, 1U);
}

static void test_coefficients_and_dwt_boundary(void) {
    CHECK(length_haar == 2U);
    check_close(haar[0], 0.707106781, "Haar coefficient 0");
    check_close(haar[1], 0.707106781, "Haar coefficient 1");
    CHECK(length_coiflet_12 == 12U);
    CHECK(length_battle_lemarie_degree_30 == 30U);

    double ***cube = create_3d_zero_array(1U, 1U, 1U);
    CHECK(cube != NULL);
    if (cube != NULL) {
        const double filter[] = {1.0};
        errno = 0;
        CHECK(discrete_wavelet_transform(
            cube,
            filter,
            1U,
            1U,
            1U,
            1U,
            1U
        ) == NULL);
        CHECK(errno == ENOTSUP);
    }
    free_3d_array(cube, 1U, 1U);

    errno = 0;
    CHECK(discrete_wavelet_transform(NULL, NULL, 0U, 0U, 0U, 0U, 0U) == NULL);
    CHECK(errno == EINVAL);
}

int main(void) {
    test_allocation_and_copy_helpers();
    test_3d_helpers();
    test_circular_shifts();
    test_filtering_and_convolution();
    test_mirror_and_sampling();
    test_coefficients_and_dwt_boundary();

    if (failures != 0U) {
        fprintf(
            stderr,
            "%zu of %zu checks failed.\n",
            failures,
            checks_run
        );
        return EXIT_FAILURE;
    }
    printf("All %zu checks passed.\n", checks_run);
    return EXIT_SUCCESS;
}
