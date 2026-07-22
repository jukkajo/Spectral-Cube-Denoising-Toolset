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

static void check_haar_case(
    const double *input,
    size_t input_length,
    const double *expected_approximation,
    const double *expected_detail,
    size_t expected_coefficient_length,
    const char *label
) {
    struct HaarDwtResult *result = haar_dwt_forward(input, input_length);
    CHECK(result != NULL);
    if (result == NULL) {
        return;
    }

    CHECK(result->original_length == input_length);
    CHECK(result->coefficient_length == expected_coefficient_length);
    CHECK(result->approximation != NULL);
    CHECK(result->detail != NULL);
    CHECK(result->approximation != result->detail);
    CHECK(result->approximation != input);
    CHECK(result->detail != input);

    for (size_t index = 0U; index < expected_coefficient_length; ++index) {
        check_close(
            result->approximation[index],
            expected_approximation[index],
            label
        );
        check_close(result->detail[index], expected_detail[index], label);
    }

    double *reconstructed = haar_dwt_inverse(result);
    CHECK(reconstructed != NULL);
    if (reconstructed != NULL) {
        double maximum_error = 0.0;
        for (size_t index = 0U; index < input_length; ++index) {
            const double error = fabs(reconstructed[index] - input[index]);
            if (error > maximum_error) {
                maximum_error = error;
            }
        }
        CHECK(maximum_error <= 1.0e-12);
    }

    free(reconstructed);
    haar_dwt_result_free(result);
}

static void test_haar_dwt(void) {
    const double sqrt_two = sqrt(2.0);
    const double inverse_sqrt_two = 1.0 / sqrt_two;

    const double equal_pair[] = {1.0, 1.0};
    const double equal_approximation[] = {sqrt_two};
    const double equal_detail[] = {0.0};
    check_haar_case(
        equal_pair,
        2U,
        equal_approximation,
        equal_detail,
        1U,
        "Haar [1, 1]"
    );

    const double opposite_pair[] = {1.0, -1.0};
    const double opposite_approximation[] = {0.0};
    const double opposite_detail[] = {sqrt_two};
    check_haar_case(
        opposite_pair,
        2U,
        opposite_approximation,
        opposite_detail,
        1U,
        "Haar [1, -1]"
    );

    const double four_values[] = {1.0, 2.0, 3.0, 4.0};
    const double four_approximation[] = {
        3.0 * inverse_sqrt_two,
        7.0 * inverse_sqrt_two
    };
    const double four_detail[] = {
        -inverse_sqrt_two,
        -inverse_sqrt_two
    };
    check_haar_case(
        four_values,
        4U,
        four_approximation,
        four_detail,
        2U,
        "Haar [1, 2, 3, 4]"
    );

    const double constant[] = {3.0, 3.0, 3.0, 3.0};
    const double constant_approximation[] = {
        3.0 * sqrt_two,
        3.0 * sqrt_two
    };
    const double constant_detail[] = {0.0, 0.0};
    check_haar_case(
        constant,
        4U,
        constant_approximation,
        constant_detail,
        2U,
        "Haar constant"
    );

    const double impulse[] = {1.0, 0.0, 0.0, 0.0};
    const double impulse_approximation[] = {inverse_sqrt_two, 0.0};
    const double impulse_detail[] = {inverse_sqrt_two, 0.0};
    check_haar_case(
        impulse,
        4U,
        impulse_approximation,
        impulse_detail,
        2U,
        "Haar impulse"
    );

    const double alternating[] = {1.0, -1.0, 1.0, -1.0};
    const double alternating_approximation[] = {0.0, 0.0};
    const double alternating_detail[] = {sqrt_two, sqrt_two};
    check_haar_case(
        alternating,
        4U,
        alternating_approximation,
        alternating_detail,
        2U,
        "Haar alternating signs"
    );

    const double ramp[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    const double ramp_approximation[] = {
        inverse_sqrt_two,
        5.0 * inverse_sqrt_two,
        9.0 * inverse_sqrt_two
    };
    const double ramp_detail[] = {
        -inverse_sqrt_two,
        -inverse_sqrt_two,
        -inverse_sqrt_two
    };
    check_haar_case(
        ramp,
        6U,
        ramp_approximation,
        ramp_detail,
        3U,
        "Haar ramp"
    );

    const double singleton[] = {5.0};
    const double singleton_approximation[] = {5.0 * sqrt_two};
    const double singleton_detail[] = {0.0};
    check_haar_case(
        singleton,
        1U,
        singleton_approximation,
        singleton_detail,
        1U,
        "Haar singleton"
    );

    const double odd[] = {1.0, 2.0, 3.0};
    const double odd_approximation[] = {
        3.0 * inverse_sqrt_two,
        3.0 * sqrt_two
    };
    const double odd_detail[] = {-inverse_sqrt_two, 0.0};
    check_haar_case(
        odd,
        3U,
        odd_approximation,
        odd_detail,
        2U,
        "Haar odd length"
    );

    double mutable_input[] = {2.0, 4.0, 6.0};
    const double original_input[] = {2.0, 4.0, 6.0};
    struct HaarDwtResult *owned = haar_dwt_forward(mutable_input, 3U);
    CHECK(owned != NULL);
    if (owned != NULL) {
        for (size_t index = 0U; index < 3U; ++index) {
            check_close(
                mutable_input[index],
                original_input[index],
                "Haar input immutability"
            );
        }

        const double original_detail = owned->detail[0];
        owned->approximation[0] = 99.0;
        check_close(
            owned->detail[0],
            original_detail,
            "Haar coefficient array ownership"
        );
        check_close(
            mutable_input[0],
            original_input[0],
            "Haar coefficient/input ownership"
        );

        double *independent_inverse = haar_dwt_inverse(owned);
        CHECK(independent_inverse != NULL);
        if (independent_inverse != NULL) {
            const double saved_approximation = owned->approximation[0];
            independent_inverse[0] = -100.0;
            check_close(
                owned->approximation[0],
                saved_approximation,
                "Haar inverse output ownership"
            );
        }
        free(independent_inverse);
    }
    haar_dwt_result_free(owned);

    errno = 0;
    CHECK(haar_dwt_forward(NULL, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_dwt_forward(equal_pair, 0U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_dwt_forward(equal_pair, SIZE_MAX) == NULL);
    CHECK(errno == EOVERFLOW);

    double coefficient = 0.0;
    struct HaarDwtResult invalid = {0U, 1U, &coefficient, &coefficient};
    errno = 0;
    CHECK(haar_dwt_inverse(NULL) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_dwt_inverse(&invalid) == NULL);
    CHECK(errno == EINVAL);

    invalid.original_length = 1U;
    invalid.coefficient_length = 0U;
    errno = 0;
    CHECK(haar_dwt_inverse(&invalid) == NULL);
    CHECK(errno == EINVAL);

    invalid.original_length = 4U;
    invalid.coefficient_length = 1U;
    errno = 0;
    CHECK(haar_dwt_inverse(&invalid) == NULL);
    CHECK(errno == EINVAL);

    invalid.original_length = 2U;
    invalid.coefficient_length = 1U;
    invalid.approximation = NULL;
    errno = 0;
    CHECK(haar_dwt_inverse(&invalid) == NULL);
    CHECK(errno == EINVAL);

    invalid.approximation = &coefficient;
    invalid.detail = NULL;
    errno = 0;
    CHECK(haar_dwt_inverse(&invalid) == NULL);
    CHECK(errno == EINVAL);

    invalid.original_length = SIZE_MAX;
    invalid.coefficient_length = SIZE_MAX / 2U + SIZE_MAX % 2U;
    invalid.detail = &coefficient;
    errno = 0;
    CHECK(haar_dwt_inverse(&invalid) == NULL);
    CHECK(errno == EOVERFLOW);

    haar_dwt_result_free(NULL);
}

static size_t test_half_rounded_up(size_t length) {
    return length / 2U + length % 2U;
}

static void fill_multilevel_pattern(
    double *values,
    size_t length,
    size_t pattern
) {
    static const double pseudo_random_values[] = {
        0.25, -1.5, 3.75, 2.125, -4.5, 0.875, 6.25, -2.75,
        1.625, -0.375, 5.5, -3.125, 2.875, 0.5, -6.75, 4.25
    };

    for (size_t index = 0U; index < length; ++index) {
        switch (pattern) {
            case 0U:
                values[index] = 3.0;
                break;
            case 1U:
                values[index] = index == 0U ? 1.0 : 0.0;
                break;
            case 2U:
                values[index] = (double)index;
                break;
            case 3U:
                values[index] = index % 2U == 0U ? 1.0 : -1.0;
                break;
            default:
                values[index] = pseudo_random_values[index];
                break;
        }
    }
}

static void check_multilevel_case(
    const double *input,
    size_t input_length,
    size_t levels,
    const char *label
) {
    double input_copy[16];
    CHECK(input_length <= 16U);
    if (input_length > 16U) {
        return;
    }
    for (size_t index = 0U; index < input_length; ++index) {
        input_copy[index] = input[index];
    }

    struct HaarDwtDecomposition *decomposition =
        haar_dwt_multilevel_forward(input, input_length, levels);
    CHECK(decomposition != NULL);
    if (decomposition == NULL) {
        return;
    }

    CHECK(decomposition->original_length == input_length);
    CHECK(decomposition->levels == levels);
    CHECK(decomposition->final_approximation != NULL);
    CHECK(decomposition->details != NULL);
    CHECK(decomposition->detail_lengths != NULL);
    CHECK(decomposition->level_input_lengths != NULL);

    size_t expected_input_length = input_length;
    for (size_t level = 0U; level < levels; ++level) {
        const size_t expected_detail_length =
            test_half_rounded_up(expected_input_length);
        CHECK(decomposition->level_input_lengths[level] ==
            expected_input_length);
        CHECK(decomposition->detail_lengths[level] ==
            expected_detail_length);
        CHECK(decomposition->details[level] != NULL);
        CHECK(decomposition->details[level] != input);
        CHECK(decomposition->details[level] !=
            decomposition->final_approximation);
        for (size_t earlier_level = 0U;
             earlier_level < level;
             ++earlier_level) {
            CHECK(decomposition->details[level] !=
                decomposition->details[earlier_level]);
        }
        expected_input_length = expected_detail_length;
    }
    CHECK(decomposition->final_approximation_length ==
        expected_input_length);
    CHECK(decomposition->final_approximation != input);

    for (size_t index = 0U; index < input_length; ++index) {
        check_close(input[index], input_copy[index], label);
    }

    if (levels == 1U) {
        struct HaarDwtResult *one_level =
            haar_dwt_forward(input, input_length);
        CHECK(one_level != NULL);
        if (one_level != NULL) {
            CHECK(one_level->coefficient_length ==
                decomposition->final_approximation_length);
            for (size_t index = 0U;
                 index < one_level->coefficient_length;
                 ++index) {
                check_close(
                    decomposition->final_approximation[index],
                    one_level->approximation[index],
                    "multilevel level-one approximation"
                );
                check_close(
                    decomposition->details[0][index],
                    one_level->detail[index],
                    "multilevel level-one detail"
                );
            }
        }
        haar_dwt_result_free(one_level);
    }

    double detail_first_values[4];
    CHECK(levels <= 4U);
    for (size_t level = 0U; level < levels && level < 4U; ++level) {
        detail_first_values[level] = decomposition->details[level][0];
    }
    const double final_first_value = decomposition->final_approximation[0];

    double *reconstructed = haar_dwt_multilevel_inverse(decomposition);
    CHECK(reconstructed != NULL);
    if (reconstructed != NULL) {
        CHECK(reconstructed != input);
        CHECK(reconstructed != decomposition->final_approximation);
        double maximum_error = 0.0;
        for (size_t index = 0U; index < input_length; ++index) {
            const double error = fabs(reconstructed[index] - input[index]);
            if (error > maximum_error) {
                maximum_error = error;
            }
        }
        CHECK(maximum_error <= 1.0e-12);

        reconstructed[0] = -999.0;
        check_close(
            decomposition->final_approximation[0],
            final_first_value,
            "multilevel inverse output ownership"
        );
    }

    check_close(
        decomposition->final_approximation[0],
        final_first_value,
        "multilevel inverse approximation immutability"
    );
    for (size_t level = 0U; level < levels && level < 4U; ++level) {
        check_close(
            decomposition->details[level][0],
            detail_first_values[level],
            "multilevel inverse detail immutability"
        );
    }

    free(reconstructed);
    haar_dwt_decomposition_free(decomposition);
}

static void test_multilevel_haar_dwt(void) {
    const size_t input_lengths[] = {1U, 2U, 3U, 4U, 7U, 8U, 15U, 16U};
    const size_t expected_maximum_levels[] = {
        0U, 1U, 2U, 2U, 3U, 3U, 4U, 4U
    };

    CHECK(haar_dwt_max_levels(0U) == 0U);
    for (size_t length_index = 0U;
         length_index < sizeof(input_lengths) / sizeof(input_lengths[0]);
         ++length_index) {
        const size_t input_length = input_lengths[length_index];
        const size_t maximum_levels =
            expected_maximum_levels[length_index];
        CHECK(haar_dwt_max_levels(input_length) == maximum_levels);

        double input[16];
        fill_multilevel_pattern(input, input_length, 4U);

        errno = 0;
        CHECK(haar_dwt_multilevel_forward(input, input_length, 0U) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(haar_dwt_multilevel_forward(
            input,
            input_length,
            maximum_levels + 1U
        ) == NULL);
        CHECK(errno == EINVAL);

        for (size_t pattern = 0U; pattern < 5U; ++pattern) {
            fill_multilevel_pattern(input, input_length, pattern);
            for (size_t levels = 1U;
                 levels <= maximum_levels;
                 ++levels) {
                check_multilevel_case(
                    input,
                    input_length,
                    levels,
                    "multilevel input immutability"
                );
            }
        }
    }

    const double pair[] = {1.0, 2.0};
    errno = 0;
    CHECK(haar_dwt_multilevel_forward(NULL, 2U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_dwt_multilevel_forward(pair, 0U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_dwt_multilevel_forward(pair, 2U, SIZE_MAX) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_dwt_multilevel_forward(pair, SIZE_MAX, 1U) == NULL);
    CHECK(errno == EOVERFLOW);

    errno = 0;
    CHECK(haar_dwt_multilevel_inverse(NULL) == NULL);
    CHECK(errno == EINVAL);

    const double metadata_input[] = {
        1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0
    };
    struct HaarDwtDecomposition *metadata =
        haar_dwt_multilevel_forward(metadata_input, 7U, 3U);
    CHECK(metadata != NULL);
    if (metadata != NULL) {
        const size_t saved_levels = metadata->levels;
        metadata->levels = 0U;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->levels = saved_levels;

        metadata->levels = 4U;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->levels = saved_levels;

        const size_t saved_original_length = metadata->original_length;
        metadata->original_length = 1U;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->original_length = saved_original_length;

        double *saved_final_approximation =
            metadata->final_approximation;
        metadata->final_approximation = NULL;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->final_approximation = saved_final_approximation;

        const size_t saved_final_length =
            metadata->final_approximation_length;
        metadata->final_approximation_length = saved_final_length + 1U;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->final_approximation_length = saved_final_length;

        double **saved_details = metadata->details;
        metadata->details = NULL;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->details = saved_details;

        size_t *saved_detail_lengths = metadata->detail_lengths;
        metadata->detail_lengths = NULL;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->detail_lengths = saved_detail_lengths;

        size_t *saved_input_lengths = metadata->level_input_lengths;
        metadata->level_input_lengths = NULL;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->level_input_lengths = saved_input_lengths;

        const size_t saved_level_input_length =
            metadata->level_input_lengths[1];
        metadata->level_input_lengths[1] = saved_level_input_length + 1U;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->level_input_lengths[1] = saved_level_input_length;

        const size_t saved_detail_length = metadata->detail_lengths[1];
        metadata->detail_lengths[1] = saved_detail_length + 1U;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->detail_lengths[1] = saved_detail_length;

        double *saved_detail = metadata->details[1];
        metadata->details[1] = NULL;
        errno = 0;
        CHECK(haar_dwt_multilevel_inverse(metadata) == NULL);
        CHECK(errno == EINVAL);
        metadata->details[1] = saved_detail;
    }
    haar_dwt_decomposition_free(metadata);

    double coefficient = 0.0;
    double *huge_details[] = {&coefficient};
    const size_t huge_coefficient_length =
        SIZE_MAX / 2U + SIZE_MAX % 2U;
    size_t huge_detail_lengths[] = {huge_coefficient_length};
    size_t huge_input_lengths[] = {SIZE_MAX};
    struct HaarDwtDecomposition huge = {
        SIZE_MAX,
        1U,
        &coefficient,
        huge_coefficient_length,
        huge_details,
        huge_detail_lengths,
        huge_input_lengths
    };
    errno = 0;
    CHECK(haar_dwt_multilevel_inverse(&huge) == NULL);
    CHECK(errno == EOVERFLOW);

    struct HaarDwtDecomposition *partial = calloc(1U, sizeof(*partial));
    CHECK(partial != NULL);
    if (partial != NULL) {
        partial->levels = 3U;
        partial->final_approximation = create_1d_array(1U);
        partial->details = calloc(3U, sizeof(*partial->details));
        partial->detail_lengths =
            calloc(3U, sizeof(*partial->detail_lengths));
        CHECK(partial->final_approximation != NULL);
        CHECK(partial->details != NULL);
        CHECK(partial->detail_lengths != NULL);
        if (partial->details != NULL) {
            partial->details[0] = create_1d_array(1U);
            CHECK(partial->details[0] != NULL);
        }
    }
    haar_dwt_decomposition_free(partial);
    haar_dwt_decomposition_free(NULL);
}

static void test_coefficients_and_generic_dwt_boundary(void) {
    const double inverse_sqrt_two = 1.0 / sqrt(2.0);
    CHECK(length_haar == 2U);
    check_close(haar[0], inverse_sqrt_two, "Haar coefficient 0");
    check_close(haar[1], inverse_sqrt_two, "Haar coefficient 1");
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
    test_haar_dwt();
    test_multilevel_haar_dwt();
    test_coefficients_and_generic_dwt_boundary();

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
