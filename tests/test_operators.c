#define _POSIX_C_SOURCE 200809L

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "Circular-Left-Shift/circular_left_shift.h"
#include "Circular-Right-Shift/circular_right_shift.h"
#include "Command-Line/spectral_denoise_cli.h"
#include "Discrete-Wavelet-Transform/discrete_wavelet_transform.h"
#include "Haar-Denoising/haar_denoising.h"
#include "High-Pass-Downsampling-Operator/hi_pass_downsampling_operator.h"
#include "High-Pass-Upsampling-Operator/high_pass_upsampling_operator.h"
#include "Low-Pass-Downsampling-Operator/low_pass_downsampling_operator.h"
#include "Low-pass-Upsampling-Operator/low_pass_upsampling_operator.h"
#include "Mirror-Filter/mirror_filter.h"
#include "Periodic-Convolution/periodic_convolution.h"
#include "Rational-Transfer-Function/rational_transfer_function.h"
#include "Spectral-Cube-IO/spectral_cube_io.h"
#include "Spectral-Cube/spectral_cube.h"
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

static void check_vector_values(
    const double *actual,
    const double *expected,
    size_t length,
    const char *label
) {
    for (size_t index = 0U; index < length; ++index) {
        check_close(actual[index], expected[index], label);
    }
}

static void test_direct_haar_thresholding(void) {
    const double values[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    const double expected_hard[] = {-2.0, 0.0, 0.0, 0.0, 2.0};
    const double expected_soft[] = {-1.0, 0.0, 0.0, 0.0, 1.0};

    double hard_values[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    CHECK(haar_threshold_coefficients(
        hard_values,
        5U,
        1.0,
        HAAR_THRESHOLD_HARD
    ) == 0);
    check_vector_values(
        hard_values,
        expected_hard,
        5U,
        "direct hard threshold"
    );

    double soft_values[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    CHECK(haar_threshold_coefficients(
        soft_values,
        5U,
        1.0,
        HAAR_THRESHOLD_SOFT
    ) == 0);
    check_vector_values(
        soft_values,
        expected_soft,
        5U,
        "direct soft threshold"
    );

    double hard_equal[] = {-1.0, 1.0};
    const double hard_equal_expected[] = {-1.0, 1.0};
    CHECK(haar_threshold_coefficients(
        hard_equal,
        2U,
        1.0,
        HAAR_THRESHOLD_HARD
    ) == 0);
    check_vector_values(
        hard_equal,
        hard_equal_expected,
        2U,
        "hard threshold equality"
    );

    double soft_equal[] = {-1.0, 1.0};
    const double soft_equal_expected[] = {0.0, 0.0};
    CHECK(haar_threshold_coefficients(
        soft_equal,
        2U,
        1.0,
        HAAR_THRESHOLD_SOFT
    ) == 0);
    check_vector_values(
        soft_equal,
        soft_equal_expected,
        2U,
        "soft threshold equality"
    );

    double hard_zero[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    CHECK(haar_threshold_coefficients(
        hard_zero,
        5U,
        0.0,
        HAAR_THRESHOLD_HARD
    ) == 0);
    check_vector_values(hard_zero, values, 5U, "zero hard threshold");

    double soft_zero[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    CHECK(haar_threshold_coefficients(
        soft_zero,
        5U,
        0.0,
        HAAR_THRESHOLD_SOFT
    ) == 0);
    check_vector_values(soft_zero, values, 5U, "zero soft threshold");

    double hard_large[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    const double all_zero[] = {0.0, 0.0, 0.0, 0.0, 0.0};
    CHECK(haar_threshold_coefficients(
        hard_large,
        5U,
        3.0,
        HAAR_THRESHOLD_HARD
    ) == 0);
    check_vector_values(hard_large, all_zero, 5U, "large hard threshold");

    double soft_large[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    CHECK(haar_threshold_coefficients(
        soft_large,
        5U,
        3.0,
        HAAR_THRESHOLD_SOFT
    ) == 0);
    check_vector_values(soft_large, all_zero, 5U, "large soft threshold");

    double invalid_values[] = {-2.0, -0.5, 0.0, 0.5, 2.0};
    errno = 0;
    CHECK(haar_threshold_coefficients(
        invalid_values,
        5U,
        -1.0,
        HAAR_THRESHOLD_HARD
    ) == -1);
    CHECK(errno == EINVAL);
    check_vector_values(
        invalid_values,
        values,
        5U,
        "negative threshold immutability"
    );

    errno = 0;
    CHECK(haar_threshold_coefficients(
        invalid_values,
        5U,
        NAN,
        HAAR_THRESHOLD_HARD
    ) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_threshold_coefficients(
        invalid_values,
        5U,
        INFINITY,
        HAAR_THRESHOLD_SOFT
    ) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_threshold_coefficients(
        invalid_values,
        5U,
        1.0,
        (enum HaarThresholdMode)99
    ) == -1);
    CHECK(errno == EINVAL);
    check_vector_values(
        invalid_values,
        values,
        5U,
        "invalid threshold argument immutability"
    );

    errno = 0;
    CHECK(haar_threshold_coefficients(
        NULL,
        5U,
        1.0,
        HAAR_THRESHOLD_HARD
    ) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_threshold_coefficients(
        invalid_values,
        0U,
        1.0,
        HAAR_THRESHOLD_HARD
    ) == -1);
    CHECK(errno == EINVAL);
}

static void check_detail_level_selection(
    const double *input,
    size_t input_length,
    size_t levels,
    const struct HaarDetailLevelRange *range
) {
    struct HaarDwtDecomposition *baseline =
        haar_dwt_multilevel_forward(input, input_length, levels);
    struct HaarDwtDecomposition *thresholded =
        haar_dwt_multilevel_forward(input, input_length, levels);
    CHECK(baseline != NULL);
    CHECK(thresholded != NULL);
    if (baseline == NULL || thresholded == NULL) {
        haar_dwt_decomposition_free(baseline);
        haar_dwt_decomposition_free(thresholded);
        return;
    }

    CHECK(haar_threshold_detail_levels(
        thresholded,
        1000.0,
        HAAR_THRESHOLD_HARD,
        range
    ) == 0);

    check_vector_values(
        thresholded->final_approximation,
        baseline->final_approximation,
        baseline->final_approximation_length,
        "threshold approximation immutability"
    );

    for (size_t level = 0U; level < levels; ++level) {
        const int selected = range == NULL ||
            (level >= range->first_level && level <= range->last_level);
        for (size_t index = 0U;
             index < baseline->detail_lengths[level];
             ++index) {
            const double expected = selected != 0
                ? 0.0
                : baseline->details[level][index];
            check_close(
                thresholded->details[level][index],
                expected,
                "selected detail-level threshold"
            );
        }
    }

    haar_dwt_decomposition_free(baseline);
    haar_dwt_decomposition_free(thresholded);
}

static void test_haar_detail_thresholding(void) {
    const double input[] = {0.0, 1.0, 4.0, 2.0, 8.0, 3.0, -1.0, 5.0};
    const struct HaarDetailLevelRange first_level = {0U, 0U};
    const struct HaarDetailLevelRange deepest_level = {2U, 2U};
    const struct HaarDetailLevelRange contiguous_levels = {0U, 1U};

    check_detail_level_selection(input, 8U, 3U, NULL);
    check_detail_level_selection(input, 8U, 3U, &first_level);
    check_detail_level_selection(input, 8U, 3U, &deepest_level);
    check_detail_level_selection(input, 8U, 3U, &contiguous_levels);

    const double one_level_input[] = {1.0, 3.0, 2.0};
    const struct HaarDetailLevelRange only_level = {0U, 0U};
    check_detail_level_selection(
        one_level_input,
        3U,
        1U,
        &only_level
    );

    struct HaarDwtDecomposition *baseline =
        haar_dwt_multilevel_forward(input, 8U, 3U);
    struct HaarDwtDecomposition *invalid =
        haar_dwt_multilevel_forward(input, 8U, 3U);
    CHECK(baseline != NULL);
    CHECK(invalid != NULL);
    if (baseline != NULL && invalid != NULL) {
        const struct HaarDetailLevelRange reversed = {2U, 1U};
        const struct HaarDetailLevelRange out_of_range = {0U, 3U};

        errno = 0;
        CHECK(haar_threshold_detail_levels(
            invalid,
            1.0,
            HAAR_THRESHOLD_HARD,
            &reversed
        ) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(haar_threshold_detail_levels(
            invalid,
            1.0,
            HAAR_THRESHOLD_HARD,
            &out_of_range
        ) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(haar_threshold_detail_levels(
            invalid,
            -1.0,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(haar_threshold_detail_levels(
            invalid,
            NAN,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(haar_threshold_detail_levels(
            invalid,
            INFINITY,
            HAAR_THRESHOLD_SOFT,
            NULL
        ) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(haar_threshold_detail_levels(
            invalid,
            1.0,
            (enum HaarThresholdMode)99,
            NULL
        ) == -1);
        CHECK(errno == EINVAL);

        check_vector_values(
            invalid->final_approximation,
            baseline->final_approximation,
            baseline->final_approximation_length,
            "invalid detail threshold approximation immutability"
        );
        for (size_t level = 0U; level < baseline->levels; ++level) {
            check_vector_values(
                invalid->details[level],
                baseline->details[level],
                baseline->detail_lengths[level],
                "invalid detail threshold immutability"
            );
        }
    }
    haar_dwt_decomposition_free(baseline);
    haar_dwt_decomposition_free(invalid);

    errno = 0;
    CHECK(haar_threshold_detail_levels(
        NULL,
        1.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == -1);
    CHECK(errno == EINVAL);
}

static void check_zero_threshold_denoising(
    const double *input,
    size_t input_length,
    size_t levels,
    enum HaarThresholdMode mode
) {
    double input_copy[8];
    CHECK(input_length <= 8U);
    if (input_length > 8U) {
        return;
    }
    for (size_t index = 0U; index < input_length; ++index) {
        input_copy[index] = input[index];
    }

    double *output = haar_denoise_1d(
        input,
        input_length,
        levels,
        0.0,
        mode,
        NULL
    );
    CHECK(output != NULL);
    if (output != NULL) {
        CHECK(output != input);
        double maximum_error = 0.0;
        for (size_t index = 0U; index < input_length; ++index) {
            const double error = fabs(output[index] - input[index]);
            if (error > maximum_error) {
                maximum_error = error;
            }
        }
        CHECK(maximum_error <= 1.0e-12);
        output[0] = -500.0;
    }

    check_vector_values(
        input,
        input_copy,
        input_length,
        "denoising input immutability"
    );
    free(output);
}

static void test_haar_denoising(void) {
    const double one_level_odd[] = {1.0, -2.0, 4.0};
    const double multilevel_odd[] = {1.0, 3.0, -2.0, 4.0, 0.5, 6.0, -1.0};
    check_zero_threshold_denoising(
        one_level_odd,
        3U,
        1U,
        HAAR_THRESHOLD_HARD
    );
    check_zero_threshold_denoising(
        multilevel_odd,
        7U,
        3U,
        HAAR_THRESHOLD_HARD
    );
    check_zero_threshold_denoising(
        multilevel_odd,
        7U,
        3U,
        HAAR_THRESHOLD_SOFT
    );

    const double noisy_fixture[] = {1.0, 3.0, 1.0, 3.0};
    const double inverse_sqrt_two = 1.0 / sqrt(2.0);
    const double expected_soft[] = {
        1.0 + inverse_sqrt_two,
        3.0 - inverse_sqrt_two,
        1.0 + inverse_sqrt_two,
        3.0 - inverse_sqrt_two
    };
    double *soft_output = haar_denoise_1d(
        noisy_fixture,
        4U,
        2U,
        1.0,
        HAAR_THRESHOLD_SOFT,
        NULL
    );
    CHECK(soft_output != NULL);
    if (soft_output != NULL) {
        check_vector_values(
            soft_output,
            expected_soft,
            4U,
            "deterministic soft denoising fixture"
        );
    }
    free(soft_output);

    const double range_fixture[] = {1.0, 3.0, 2.0, 4.0};
    const struct HaarDetailLevelRange first_level = {0U, 0U};
    const struct HaarDetailLevelRange deepest_level = {1U, 1U};
    const double expected_all[] = {2.5, 2.5, 2.5, 2.5};
    const double expected_first[] = {2.0, 2.0, 3.0, 3.0};
    const double expected_deepest[] = {1.5, 3.5, 1.5, 3.5};

    double *all_output = haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        100.0,
        HAAR_THRESHOLD_HARD,
        NULL
    );
    double *first_output = haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        100.0,
        HAAR_THRESHOLD_HARD,
        &first_level
    );
    double *deepest_output = haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        100.0,
        HAAR_THRESHOLD_HARD,
        &deepest_level
    );
    CHECK(all_output != NULL);
    CHECK(first_output != NULL);
    CHECK(deepest_output != NULL);
    if (all_output != NULL) {
        check_vector_values(
            all_output,
            expected_all,
            4U,
            "all-level denoising"
        );
    }
    if (first_output != NULL) {
        check_vector_values(
            first_output,
            expected_first,
            4U,
            "first-level denoising"
        );
    }
    if (deepest_output != NULL) {
        check_vector_values(
            deepest_output,
            expected_deepest,
            4U,
            "deepest-level denoising"
        );
    }
    free(all_output);
    free(first_output);
    free(deepest_output);

    const double singleton[] = {5.0};
    errno = 0;
    CHECK(haar_denoise_1d(
        singleton,
        1U,
        1U,
        0.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        singleton,
        1U,
        0U,
        0.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);

    const struct HaarDetailLevelRange reversed = {1U, 0U};
    const struct HaarDetailLevelRange out_of_range = {0U, 2U};
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        1.0,
        HAAR_THRESHOLD_HARD,
        &reversed
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        1.0,
        HAAR_THRESHOLD_HARD,
        &out_of_range
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        -1.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        NAN,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        INFINITY,
        HAAR_THRESHOLD_SOFT,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        2U,
        1.0,
        (enum HaarThresholdMode)99,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        NULL,
        4U,
        2U,
        1.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        0U,
        1U,
        1.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        4U,
        3U,
        1.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(haar_denoise_1d(
        range_fixture,
        SIZE_MAX,
        1U,
        0.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EOVERFLOW);

    check_vector_values(
        noisy_fixture,
        (const double[]){1.0, 3.0, 1.0, 3.0},
        4U,
        "noisy fixture input immutability"
    );
    check_vector_values(
        range_fixture,
        (const double[]){1.0, 3.0, 2.0, 4.0},
        4U,
        "range fixture input immutability"
    );
}

static SpectralCube *cube_from_values(
    size_t height,
    size_t width,
    size_t bands,
    const double *values
) {
    SpectralCube *cube = spectral_cube_create(height, width, bands);
    if (cube == NULL) {
        return NULL;
    }

    size_t element_count = 0U;
    if (spectral_cube_element_count(cube, &element_count) != 0) {
        spectral_cube_free(cube);
        return NULL;
    }
    for (size_t index = 0U; index < element_count; ++index) {
        cube->data[index] = values[index];
    }
    return cube;
}

static void check_cube_values(
    const SpectralCube *cube,
    const double *expected,
    const char *label
) {
    size_t element_count = 0U;
    const int count_status =
        spectral_cube_element_count(cube, &element_count);
    CHECK(count_status == 0);
    if (count_status != 0) {
        return;
    }
    check_vector_values(cube->data, expected, element_count, label);
}

static void test_spectral_cube_model(void) {
    struct CubeDimensions {
        size_t height;
        size_t width;
        size_t bands;
    };
    const struct CubeDimensions dimensions[] = {
        {1U, 1U, 1U},
        {1U, 2U, 3U},
        {2U, 2U, 4U},
        {3U, 2U, 5U}
    };

    for (size_t dimension_index = 0U;
         dimension_index < sizeof(dimensions) / sizeof(dimensions[0]);
         ++dimension_index) {
        const size_t height = dimensions[dimension_index].height;
        const size_t width = dimensions[dimension_index].width;
        const size_t bands = dimensions[dimension_index].bands;
        const size_t expected_count = height * width * bands;
        SpectralCube *cube = spectral_cube_create(height, width, bands);
        CHECK(cube != NULL);
        if (cube == NULL) {
            continue;
        }

        CHECK(cube->height == height);
        CHECK(cube->width == width);
        CHECK(cube->bands == bands);
        CHECK(cube->data != NULL);

        size_t element_count = 0U;
        CHECK(spectral_cube_element_count(cube, &element_count) == 0);
        CHECK(element_count == expected_count);
        for (size_t index = 0U; index < element_count; ++index) {
            check_close(cube->data[index], 0.0, "zero-initialized cube");
        }

        for (size_t row = 0U; row < height; ++row) {
            for (size_t column = 0U; column < width; ++column) {
                for (size_t band = 0U; band < bands; ++band) {
                    const size_t expected_index =
                        ((row * width) + column) * bands + band;
                    const double expected_value =
                        (double)(expected_index + 1U);
                    CHECK(spectral_cube_set(
                        cube,
                        row,
                        column,
                        band,
                        expected_value
                    ) == 0);
                    check_close(
                        cube->data[expected_index],
                        expected_value,
                        "pixel-major cube mapping"
                    );

                    double actual_value = 0.0;
                    CHECK(spectral_cube_get(
                        cube,
                        row,
                        column,
                        band,
                        &actual_value
                    ) == 0);
                    check_close(
                        actual_value,
                        expected_value,
                        "checked cube access"
                    );
                }
            }
        }

        double first_value = 0.0;
        double last_value = 0.0;
        CHECK(spectral_cube_get(cube, 0U, 0U, 0U, &first_value) == 0);
        CHECK(spectral_cube_get(
            cube,
            height - 1U,
            width - 1U,
            bands - 1U,
            &last_value
        ) == 0);
        check_close(first_value, 1.0, "first cube element");
        check_close(last_value, (double)element_count, "last cube element");

        SpectralCube *copy = spectral_cube_copy(cube);
        CHECK(copy != NULL);
        if (copy != NULL) {
            CHECK(copy != cube);
            CHECK(copy->data != cube->data);
            CHECK(copy->height == cube->height);
            CHECK(copy->width == cube->width);
            CHECK(copy->bands == cube->bands);
            check_vector_values(
                copy->data,
                cube->data,
                element_count,
                "spectral cube copy"
            );
            CHECK(spectral_cube_set(copy, 0U, 0U, 0U, -50.0) == 0);
            check_close(cube->data[0], 1.0, "spectral cube copy ownership");
        }
        spectral_cube_free(copy);
        spectral_cube_free(cube);
    }

    SpectralCube *bounds = spectral_cube_create(2U, 2U, 4U);
    CHECK(bounds != NULL);
    if (bounds != NULL) {
        double value = 0.0;
        errno = 0;
        CHECK(spectral_cube_get(bounds, 2U, 0U, 0U, &value) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_get(bounds, 0U, 2U, 0U, &value) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_get(bounds, 0U, 0U, 4U, &value) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_set(bounds, 2U, 0U, 0U, 1.0) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_set(bounds, 0U, 2U, 0U, 1.0) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_set(bounds, 0U, 0U, 4U, 1.0) == -1);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_get(bounds, 0U, 0U, 0U, NULL) == -1);
        CHECK(errno == EINVAL);
    }
    spectral_cube_free(bounds);

    errno = 0;
    CHECK(spectral_cube_create(0U, 1U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_create(1U, 0U, 1U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_create(1U, 1U, 0U) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_create(SIZE_MAX, 2U, 1U) == NULL);
    CHECK(errno == EOVERFLOW);
    errno = 0;
    CHECK(spectral_cube_create(1U, SIZE_MAX, 2U) == NULL);
    CHECK(errno == EOVERFLOW);
    errno = 0;
    CHECK(spectral_cube_create(
        1U,
        1U,
        SIZE_MAX / sizeof(double) + 1U
    ) == NULL);
    CHECK(errno == EOVERFLOW);

    size_t element_count = 0U;
    errno = 0;
    CHECK(spectral_cube_element_count(NULL, &element_count) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    SpectralCube invalid_cube = {1U, 1U, 1U, NULL};
    CHECK(spectral_cube_element_count(&invalid_cube, &element_count) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_copy(NULL) == NULL);
    CHECK(errno == EINVAL);

    spectral_cube_free(NULL);
}

static void check_zero_threshold_cube(
    SpectralCube *source,
    size_t levels,
    enum HaarThresholdMode mode
) {
    SpectralCube *snapshot = spectral_cube_copy(source);
    CHECK(snapshot != NULL);
    if (snapshot == NULL) {
        return;
    }

    SpectralCube *destination = spectral_cube_denoise_spectral(
        source,
        levels,
        0.0,
        mode,
        NULL
    );
    CHECK(destination != NULL);
    if (destination != NULL) {
        CHECK(destination != source);
        CHECK(destination->data != source->data);
        CHECK(destination->height == source->height);
        CHECK(destination->width == source->width);
        CHECK(destination->bands == source->bands);

        size_t element_count = 0U;
        CHECK(spectral_cube_element_count(source, &element_count) == 0);
        double maximum_error = 0.0;
        for (size_t index = 0U; index < element_count; ++index) {
            const double error =
                fabs(destination->data[index] - source->data[index]);
            if (error > maximum_error) {
                maximum_error = error;
            }
            CHECK(isfinite(destination->data[index]));
        }
        CHECK(maximum_error <= 1.0e-12);

        destination->data[0] = -1000.0;
        check_close(
            source->data[0],
            snapshot->data[0],
            "spectral denoising destination ownership"
        );
    }

    size_t source_count = 0U;
    CHECK(spectral_cube_element_count(source, &source_count) == 0);
    check_vector_values(
        source->data,
        snapshot->data,
        source_count,
        "spectral denoising source immutability"
    );
    spectral_cube_free(destination);
    spectral_cube_free(snapshot);
}

static void check_cube_range_against_1d(
    const SpectralCube *source,
    size_t levels,
    const struct HaarDetailLevelRange *range
) {
    SpectralCube *destination = spectral_cube_denoise_spectral(
        source,
        levels,
        100.0,
        HAAR_THRESHOLD_HARD,
        range
    );
    CHECK(destination != NULL);
    if (destination == NULL) {
        return;
    }

    double *expected = haar_denoise_1d(
        source->data,
        source->bands,
        levels,
        100.0,
        HAAR_THRESHOLD_HARD,
        range
    );
    CHECK(expected != NULL);
    if (expected != NULL) {
        check_vector_values(
            destination->data,
            expected,
            source->bands,
            "spectral cube level range"
        );
    }
    free(expected);
    spectral_cube_free(destination);
}

static void test_spectral_cube_denoising(void) {
    const double two_band_values[] = {1.0, 3.0};
    SpectralCube *two_band = cube_from_values(1U, 1U, 2U, two_band_values);
    CHECK(two_band != NULL);
    if (two_band != NULL) {
        SpectralCube *hard = spectral_cube_denoise_spectral(
            two_band,
            1U,
            2.0,
            HAAR_THRESHOLD_HARD,
            NULL
        );
        CHECK(hard != NULL);
        if (hard != NULL) {
            check_cube_values(
                hard,
                (const double[]){2.0, 2.0},
                "direct two-band hard denoising"
            );
        }
        spectral_cube_free(hard);

        const double inverse_sqrt_two = 1.0 / sqrt(2.0);
        SpectralCube *soft = spectral_cube_denoise_spectral(
            two_band,
            1U,
            1.0,
            HAAR_THRESHOLD_SOFT,
            NULL
        );
        CHECK(soft != NULL);
        if (soft != NULL) {
            const double expected_soft[] = {
                1.0 + inverse_sqrt_two,
                3.0 - inverse_sqrt_two
            };
            check_cube_values(
                soft,
                expected_soft,
                "direct two-band soft denoising"
            );
        }
        spectral_cube_free(soft);
    }
    spectral_cube_free(two_band);

    const double multiple_pixel_values[] = {
        1.0, 3.0, 2.0, 4.0,
        0.0, 2.0, 4.0, 6.0
    };
    const double multiple_pixel_expected[] = {
        2.5, 2.5, 2.5, 2.5,
        3.0, 3.0, 3.0, 3.0
    };
    SpectralCube *multiple_pixels = cube_from_values(
        1U,
        2U,
        4U,
        multiple_pixel_values
    );
    CHECK(multiple_pixels != NULL);
    if (multiple_pixels != NULL) {
        SpectralCube *denoised = spectral_cube_denoise_spectral(
            multiple_pixels,
            2U,
            100.0,
            HAAR_THRESHOLD_HARD,
            NULL
        );
        CHECK(denoised != NULL);
        if (denoised != NULL) {
            check_cube_values(
                denoised,
                multiple_pixel_expected,
                "independent pixel spectra"
            );
            CHECK(denoised->height == 1U);
            CHECK(denoised->width == 2U);
            CHECK(denoised->bands == 4U);
        }
        check_cube_values(
            multiple_pixels,
            multiple_pixel_values,
            "direct cube source immutability"
        );
        spectral_cube_free(denoised);
    }

    check_zero_threshold_cube(
        multiple_pixels,
        2U,
        HAAR_THRESHOLD_HARD
    );
    spectral_cube_free(multiple_pixels);

    double odd_values[30];
    for (size_t index = 0U; index < 30U; ++index) {
        odd_values[index] = (double)index * 0.25 - 2.0;
    }
    SpectralCube *odd_bands = cube_from_values(3U, 2U, 5U, odd_values);
    CHECK(odd_bands != NULL);
    if (odd_bands != NULL) {
        check_zero_threshold_cube(
            odd_bands,
            3U,
            HAAR_THRESHOLD_SOFT
        );
    }
    spectral_cube_free(odd_bands);

    const double range_values[] = {
        0.0, 1.0, 4.0, 2.0, 8.0, 3.0, -1.0, 5.0
    };
    SpectralCube *range_cube = cube_from_values(
        1U,
        1U,
        8U,
        range_values
    );
    CHECK(range_cube != NULL);
    if (range_cube != NULL) {
        const struct HaarDetailLevelRange first = {0U, 0U};
        const struct HaarDetailLevelRange deepest = {2U, 2U};
        const struct HaarDetailLevelRange contiguous = {0U, 1U};
        check_cube_range_against_1d(range_cube, 3U, &first);
        check_cube_range_against_1d(range_cube, 3U, &deepest);
        check_cube_range_against_1d(range_cube, 3U, &contiguous);
    }
    spectral_cube_free(range_cube);

    SpectralCube *invalid_source = cube_from_values(
        1U,
        2U,
        4U,
        multiple_pixel_values
    );
    SpectralCube *invalid_snapshot = spectral_cube_copy(invalid_source);
    CHECK(invalid_source != NULL);
    CHECK(invalid_snapshot != NULL);
    if (invalid_source != NULL) {
        const struct HaarDetailLevelRange reversed = {1U, 0U};
        const struct HaarDetailLevelRange out_of_range = {0U, 2U};
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            2U,
            1.0,
            HAAR_THRESHOLD_HARD,
            &reversed
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            2U,
            1.0,
            HAAR_THRESHOLD_HARD,
            &out_of_range
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            0U,
            1.0,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            3U,
            1.0,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            2U,
            -1.0,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            2U,
            NAN,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            2U,
            INFINITY,
            HAAR_THRESHOLD_SOFT,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            invalid_source,
            2U,
            1.0,
            (enum HaarThresholdMode)99,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
    }
    if (invalid_source != NULL && invalid_snapshot != NULL) {
        check_cube_values(
            invalid_source,
            invalid_snapshot->data,
            "failed denoising source immutability"
        );
    }
    spectral_cube_free(invalid_snapshot);
    spectral_cube_free(invalid_source);

    SpectralCube *single_band = spectral_cube_create(1U, 1U, 1U);
    CHECK(single_band != NULL);
    if (single_band != NULL) {
        errno = 0;
        CHECK(spectral_cube_denoise_spectral(
            single_band,
            1U,
            0.0,
            HAAR_THRESHOLD_HARD,
            NULL
        ) == NULL);
        CHECK(errno == EINVAL);
    }
    spectral_cube_free(single_band);

    errno = 0;
    CHECK(spectral_cube_denoise_spectral(
        NULL,
        1U,
        0.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EINVAL);

    double placeholder = 0.0;
    SpectralCube overflow_cube = {
        1U,
        1U,
        SIZE_MAX / sizeof(double) + 1U,
        &placeholder
    };
    errno = 0;
    CHECK(spectral_cube_denoise_spectral(
        &overflow_cube,
        1U,
        0.0,
        HAAR_THRESHOLD_HARD,
        NULL
    ) == NULL);
    CHECK(errno == EOVERFLOW);
}

static void test_encode_uint64_le(
    uint64_t value,
    unsigned char bytes[8]
) {
    for (size_t index = 0U; index < 8U; ++index) {
        bytes[index] = (unsigned char)(value & UINT64_C(0xff));
        value >>= 8U;
    }
}

static void test_encode_double_le(
    double value,
    unsigned char bytes[8]
) {
    uint64_t bits = 0U;
    memcpy(&bits, &value, sizeof(bits));
    test_encode_uint64_le(bits, bytes);
}

static int write_raw_file(
    const char *path,
    const unsigned char *bytes,
    size_t length
) {
    FILE *stream = fopen(path, "wb");
    if (stream == NULL) {
        return -1;
    }
    const size_t written = fwrite(bytes, 1U, length, stream);
    if (written != length || fclose(stream) != 0) {
        if (written != length) {
            (void)fclose(stream);
        }
        return -1;
    }
    return 0;
}

static int write_raw_cube(
    const char magic[8],
    uint64_t height,
    uint64_t width,
    uint64_t bands,
    const unsigned char *payload,
    size_t payload_length,
    const char *path
) {
    unsigned char header[SPECTRAL_CUBE_FILE_HEADER_SIZE] = {0U};
    memcpy(header, magic, 8U);
    test_encode_uint64_le(height, header + 8U);
    test_encode_uint64_le(width, header + 16U);
    test_encode_uint64_le(bands, header + 24U);

    FILE *stream = fopen(path, "wb");
    if (stream == NULL) {
        return -1;
    }
    int status = 0;
    if (fwrite(header, 1U, sizeof(header), stream) != sizeof(header)) {
        status = -1;
    } else if (payload_length != 0U &&
               fwrite(payload, 1U, payload_length, stream) != payload_length) {
        status = -1;
    }
    if (fclose(stream) != 0) {
        status = -1;
    }
    return status;
}

static int test_path_exists(const char *path) {
    FILE *stream = fopen(path, "rb");
    if (stream == NULL) {
        return 0;
    }
    (void)fclose(stream);
    return 1;
}

static int test_temporary_paths_absent(const char *destination_path) {
    char path[512];
    for (unsigned int attempt = 0U; attempt < 100U; ++attempt) {
        const int written = snprintf(
            path,
            sizeof(path),
            "%s.tmp.%u",
            destination_path,
            attempt
        );
        if (written < 0 || (size_t)written >= sizeof(path) ||
            test_path_exists(path) != 0) {
            return 0;
        }
    }
    return 1;
}

static int read_file_bytes(
    const char *path,
    unsigned char *bytes,
    size_t capacity,
    size_t *length
) {
    FILE *stream = fopen(path, "rb");
    if (stream == NULL) {
        return -1;
    }
    const size_t amount = fread(bytes, 1U, capacity, stream);
    if (ferror(stream) != 0 || fclose(stream) != 0) {
        return -1;
    }
    *length = amount;
    return 0;
}

static void check_read_failure(const char *path, int expected_errno) {
    errno = 0;
    SpectralCube *cube = spectral_cube_read_file(path);
    CHECK(cube == NULL);
    CHECK(errno == expected_errno);
    spectral_cube_free(cube);
}

static void test_spectral_cube_file_io(void) {
    const char *valid_path = "tests/.tmp_valid_cube.scube";
    const char *malformed_path = "tests/.tmp_malformed_cube.scube";
    const char *existing_path = "tests/.tmp_existing_cube.scube";
    const char *existing_symlink =
        "tests/.tmp_existing_cube_symlink.scube";
    const char *missing_symlink_target =
        "tests/.tmp_missing_symlink_target.scube";
    const char *not_directory = "tests/.tmp_not_a_directory";
    const char *invalid_output =
        "tests/.tmp_not_a_directory/output.scube";
    const double values[] = {1.0, -2.0, 3.5, 4.0, 5.25, -6.0};
    unsigned char payload[16] = {0U};
    unsigned char file_bytes[96] = {0U};
    size_t file_length = 0U;

    (void)remove(valid_path);
    (void)remove(malformed_path);
    (void)remove(existing_path);
    (void)unlink(existing_symlink);
    (void)remove(missing_symlink_target);
    (void)remove(not_directory);

    SpectralCube *source = cube_from_values(1U, 2U, 3U, values);
    CHECK(source != NULL);
    if (source == NULL) {
        return;
    }

    CHECK(spectral_cube_write_file(valid_path, source) == 0);
    CHECK(test_path_exists(valid_path) != 0);
    CHECK(test_temporary_paths_absent(valid_path) != 0);
    check_cube_values(source, values, "file write source immutability");

    SpectralCube *round_trip = spectral_cube_read_file(valid_path);
    CHECK(round_trip != NULL);
    if (round_trip != NULL) {
        CHECK(round_trip->height == 1U);
        CHECK(round_trip->width == 2U);
        CHECK(round_trip->bands == 3U);
        CHECK(round_trip->data != source->data);
        check_cube_values(round_trip, values, "cube file exact round trip");
        round_trip->data[0] = 99.0;
        check_close(source->data[0], 1.0, "cube file independent ownership");
    }
    spectral_cube_free(round_trip);

    CHECK(read_file_bytes(
        valid_path,
        file_bytes,
        sizeof(file_bytes),
        &file_length
    ) == 0);
    CHECK(file_length == 32U + sizeof(values));
    CHECK(memcmp(file_bytes, SPECTRAL_CUBE_FILE_MAGIC, 8U) == 0);
    unsigned char expected_bytes[8] = {0U};
    test_encode_uint64_le(UINT64_C(1), expected_bytes);
    CHECK(memcmp(file_bytes + 8U, expected_bytes, 8U) == 0);
    test_encode_uint64_le(UINT64_C(2), expected_bytes);
    CHECK(memcmp(file_bytes + 16U, expected_bytes, 8U) == 0);
    test_encode_uint64_le(UINT64_C(3), expected_bytes);
    CHECK(memcmp(file_bytes + 24U, expected_bytes, 8U) == 0);
    for (size_t index = 0U; index < 6U; ++index) {
        test_encode_double_le(values[index], expected_bytes);
        CHECK(memcmp(file_bytes + 32U + index * 8U, expected_bytes, 8U) == 0);
    }

    CHECK(write_raw_cube(
        "XCUBE001",
        UINT64_C(1),
        UINT64_C(1),
        UINT64_C(1),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);
    CHECK(write_raw_cube(
        "SCUBE002",
        UINT64_C(1),
        UINT64_C(1),
        UINT64_C(1),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);

    CHECK(write_raw_file(
        malformed_path,
        (const unsigned char *)SPECTRAL_CUBE_FILE_MAGIC,
        7U
    ) == 0);
    check_read_failure(malformed_path, EINVAL);

    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_C(0),
        UINT64_C(1),
        UINT64_C(1),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);
    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_C(1),
        UINT64_C(0),
        UINT64_C(1),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);
    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_C(1),
        UINT64_C(1),
        UINT64_C(0),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);

    const uint64_t conversion_dimension = UINT64_MAX;
    if ((uint64_t)(size_t)conversion_dimension != conversion_dimension) {
        CHECK(write_raw_cube(
            SPECTRAL_CUBE_FILE_MAGIC,
            conversion_dimension,
            UINT64_C(1),
            UINT64_C(1),
            NULL,
            0U,
            malformed_path
        ) == 0);
        check_read_failure(malformed_path, EOVERFLOW);
    }

    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_MAX,
        UINT64_C(2),
        UINT64_C(1),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EOVERFLOW);
    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_C(1),
        UINT64_C(1),
        (uint64_t)(SIZE_MAX / sizeof(double)) + UINT64_C(1),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EOVERFLOW);

    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_C(1),
        UINT64_C(1),
        UINT64_C(2),
        NULL,
        0U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);
    test_encode_double_le(1.0, payload);
    test_encode_double_le(2.0, payload + 8U);
    CHECK(write_raw_cube(
        SPECTRAL_CUBE_FILE_MAGIC,
        UINT64_C(1),
        UINT64_C(1),
        UINT64_C(2),
        payload,
        sizeof(payload) - 1U,
        malformed_path
    ) == 0);
    check_read_failure(malformed_path, EINVAL);

    file_bytes[file_length] = 0xffU;
    CHECK(write_raw_file(
        malformed_path,
        file_bytes,
        file_length + 1U
    ) == 0);
    check_read_failure(malformed_path, EINVAL);

    (void)remove("tests/.tmp_missing_cube.scube");
    errno = 0;
    CHECK(spectral_cube_read_file("tests/.tmp_missing_cube.scube") == NULL);
    CHECK(errno == ENOENT);

    const unsigned char marker[] = {0x41U, 0x42U, 0x43U};
    CHECK(write_raw_file(existing_path, marker, sizeof(marker)) == 0);
    errno = 0;
    CHECK(spectral_cube_write_file(existing_path, source) == -1);
    CHECK(errno == EEXIST);
    CHECK(read_file_bytes(
        existing_path,
        file_bytes,
        sizeof(file_bytes),
        &file_length
    ) == 0);
    CHECK(file_length == sizeof(marker));
    CHECK(memcmp(file_bytes, marker, sizeof(marker)) == 0);
    CHECK(test_temporary_paths_absent(existing_path) != 0);

    const char symlink_target[] = ".tmp_missing_symlink_target.scube";
    CHECK(symlink(symlink_target, existing_symlink) == 0);
    errno = 0;
    CHECK(spectral_cube_write_file(existing_symlink, source) == -1);
    CHECK(errno == EEXIST);
    char actual_symlink_target[128] = {0};
    const ssize_t target_length = readlink(
        existing_symlink,
        actual_symlink_target,
        sizeof(actual_symlink_target)
    );
    CHECK(target_length == (ssize_t)(sizeof(symlink_target) - 1U));
    if (target_length > 0 &&
        (size_t)target_length < sizeof(actual_symlink_target)) {
        actual_symlink_target[(size_t)target_length] = '\0';
        CHECK(strcmp(actual_symlink_target, symlink_target) == 0);
    }
    CHECK(test_path_exists(missing_symlink_target) == 0);
    CHECK(test_temporary_paths_absent(existing_symlink) != 0);

    CHECK(write_raw_file(not_directory, marker, sizeof(marker)) == 0);
    errno = 0;
    CHECK(spectral_cube_write_file(invalid_output, source) == -1);
    CHECK(errno == ENOTDIR);
    CHECK(test_path_exists(invalid_output) == 0);
    CHECK(test_temporary_paths_absent(invalid_output) != 0);

    SpectralCube invalid_cube = {1U, 1U, 1U, NULL};
    errno = 0;
    CHECK(spectral_cube_write_file(malformed_path, &invalid_cube) == -1);
    CHECK(errno == EINVAL);
    CHECK(test_temporary_paths_absent(malformed_path) != 0);
    double placeholder = 0.0;
    SpectralCube overflow_cube = {
        1U,
        1U,
        SIZE_MAX / sizeof(double) + 1U,
        &placeholder
    };
    errno = 0;
    CHECK(spectral_cube_write_file(malformed_path, &overflow_cube) == -1);
    CHECK(errno == EOVERFLOW);
    CHECK(test_temporary_paths_absent(malformed_path) != 0);
    errno = 0;
    CHECK(spectral_cube_write_file(NULL, source) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_write_file("", source) == -1);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_read_file(NULL) == NULL);
    CHECK(errno == EINVAL);
    errno = 0;
    CHECK(spectral_cube_read_file("") == NULL);
    CHECK(errno == EINVAL);

    check_cube_values(source, values, "file I/O source immutability");
    spectral_cube_free(source);
    (void)remove(valid_path);
    (void)remove(malformed_path);
    (void)remove(existing_path);
    (void)unlink(existing_symlink);
    (void)remove(not_directory);
}

static int run_cli_command(const char *arguments) {
    char command[2048];
    const int written = snprintf(
        command,
        sizeof(command),
        "./spectral-denoise %s > tests/.tmp_cli_stdout "
        "2> tests/.tmp_cli_stderr",
        arguments
    );
    if (written < 0 || (size_t)written >= sizeof(command)) {
        return -1;
    }
    return system(command);
}

static int file_contains_text(const char *path, const char *text) {
    unsigned char bytes[1024] = {0U};
    size_t length = 0U;
    if (read_file_bytes(path, bytes, sizeof(bytes) - 1U, &length) != 0) {
        return 0;
    }
    bytes[length] = '\0';
    return strstr((const char *)bytes, text) != NULL;
}

static void check_cli_failure(
    const char *arguments,
    const char *output_path
) {
    (void)remove(output_path);
    CHECK(run_cli_command(arguments) != 0);
    CHECK(test_path_exists(output_path) == 0);
    CHECK(test_temporary_paths_absent(output_path) != 0);
}

static void test_spectral_denoise_cli(void) {
    const char *input_path = "tests/.tmp_cli_input.scube";
    const char *hard_path = "tests/.tmp_cli_hard.scube";
    const char *soft_path = "tests/.tmp_cli_soft.scube";
    const char *failure_path = "tests/.tmp_cli_failure.scube";
    const char *existing_path = "tests/.tmp_cli_existing.scube";
    const double values[] = {
        1.0, 2.0, 4.0, 8.0,
        -1.0, 3.0, -2.0, 6.0
    };
    const unsigned char marker[] = {0x55U, 0xaaU};

    (void)remove(input_path);
    (void)remove(hard_path);
    (void)remove(soft_path);
    (void)remove(failure_path);
    (void)remove(existing_path);
    (void)remove("tests/.tmp_cli_stdout");
    (void)remove("tests/.tmp_cli_stderr");

    SpectralCube *source = cube_from_values(1U, 2U, 4U, values);
    CHECK(source != NULL);
    if (source == NULL) {
        return;
    }
    CHECK(spectral_cube_write_file(input_path, source) == 0);

    CHECK(run_cli_command("--help") == 0);
    CHECK(file_contains_text("tests/.tmp_cli_stdout", "Usage:") != 0);
    CHECK(run_cli_command("--version") == 0);
    CHECK(file_contains_text(
        "tests/.tmp_cli_stdout",
        SPECTRAL_DENOISE_VERSION
    ) != 0);

    check_cli_failure("", failure_path);
    check_cli_failure(
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube --levels 2 "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold 0.5",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold 0.5 --mode hard --unknown",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels two "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 0 "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 3 "
        "--threshold 0.5 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold -1 --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold nan --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold inf --mode hard",
        failure_path
    );
    check_cli_failure(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold 0.5 --mode garrote",
        failure_path
    );
    CHECK(run_cli_command(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_input.scube --levels 2 "
        "--threshold 0.5 --mode hard"
    ) != 0);
    CHECK(test_path_exists(input_path) != 0);
    CHECK(test_temporary_paths_absent(input_path) != 0);
    check_cli_failure(
        "--input tests/.tmp_cli_missing.scube "
        "--output tests/.tmp_cli_failure.scube --levels 2 "
        "--threshold 0.5 --mode hard",
        failure_path
    );

    CHECK(write_raw_file(existing_path, marker, sizeof(marker)) == 0);
    CHECK(run_cli_command(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_existing.scube --levels 2 "
        "--threshold 0.5 --mode hard"
    ) != 0);
    unsigned char existing_bytes[8] = {0U};
    size_t existing_length = 0U;
    CHECK(read_file_bytes(
        existing_path,
        existing_bytes,
        sizeof(existing_bytes),
        &existing_length
    ) == 0);
    CHECK(existing_length == sizeof(marker));
    CHECK(memcmp(existing_bytes, marker, sizeof(marker)) == 0);
    CHECK(test_temporary_paths_absent(existing_path) != 0);

    char *hard_arguments[] = {
        "spectral-denoise",
        "--input",
        "tests/.tmp_cli_input.scube",
        "--output",
        "tests/.tmp_cli_hard.scube",
        "--levels",
        "2",
        "--threshold",
        "0.75",
        "--mode",
        "hard"
    };
    CHECK(spectral_denoise_cli_run(11, hard_arguments) == EXIT_SUCCESS);
    SpectralCube *hard_output = spectral_cube_read_file(hard_path);
    SpectralCube *direct_hard = spectral_cube_denoise_spectral(
        source,
        2U,
        0.75,
        HAAR_THRESHOLD_HARD,
        NULL
    );
    CHECK(hard_output != NULL);
    CHECK(direct_hard != NULL);
    if (hard_output != NULL && direct_hard != NULL) {
        CHECK(hard_output->height == source->height);
        CHECK(hard_output->width == source->width);
        CHECK(hard_output->bands == source->bands);
        check_cube_values(
            hard_output,
            direct_hard->data,
            "CLI hard output matches direct API"
        );
    }
    spectral_cube_free(direct_hard);
    spectral_cube_free(hard_output);

    (void)remove(hard_path);
    CHECK(run_cli_command(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_hard.scube --levels 2 "
        "--threshold 0.75 --mode hard"
    ) == 0);
    hard_output = spectral_cube_read_file(hard_path);
    direct_hard = spectral_cube_denoise_spectral(
        source,
        2U,
        0.75,
        HAAR_THRESHOLD_HARD,
        NULL
    );
    CHECK(hard_output != NULL);
    CHECK(direct_hard != NULL);
    if (hard_output != NULL && direct_hard != NULL) {
        check_cube_values(
            hard_output,
            direct_hard->data,
            "executable hard output matches direct API"
        );
    }
    spectral_cube_free(direct_hard);
    spectral_cube_free(hard_output);

    CHECK(run_cli_command(
        "--input tests/.tmp_cli_input.scube "
        "--output tests/.tmp_cli_soft.scube --levels 2 "
        "--threshold 0.75 --mode soft"
    ) == 0);
    SpectralCube *soft_output = spectral_cube_read_file(soft_path);
    SpectralCube *direct_soft = spectral_cube_denoise_spectral(
        source,
        2U,
        0.75,
        HAAR_THRESHOLD_SOFT,
        NULL
    );
    CHECK(soft_output != NULL);
    CHECK(direct_soft != NULL);
    if (soft_output != NULL && direct_soft != NULL) {
        CHECK(soft_output->height == 1U);
        CHECK(soft_output->width == 2U);
        CHECK(soft_output->bands == 4U);
        size_t element_count = 0U;
        CHECK(spectral_cube_element_count(soft_output, &element_count) == 0);
        for (size_t index = 0U; index < element_count; ++index) {
            CHECK(isfinite(soft_output->data[index]));
        }
        check_cube_values(
            soft_output,
            direct_soft->data,
            "CLI soft output matches direct API"
        );
    }
    spectral_cube_free(direct_soft);
    spectral_cube_free(soft_output);

    SpectralCube *input_after = spectral_cube_read_file(input_path);
    CHECK(input_after != NULL);
    if (input_after != NULL) {
        check_cube_values(input_after, values, "CLI input file immutability");
    }
    spectral_cube_free(input_after);
    check_cube_values(source, values, "CLI source immutability");
    spectral_cube_free(source);

    (void)remove(input_path);
    (void)remove(hard_path);
    (void)remove(soft_path);
    (void)remove(failure_path);
    (void)remove(existing_path);
    (void)remove("tests/.tmp_cli_stdout");
    (void)remove("tests/.tmp_cli_stderr");
}

static int write_manual_cli_fixture(const char *path) {
    const double values[] = {
        1.0, 2.0, 4.0, 8.0,
        -1.0, 3.0, -2.0, 6.0
    };
    SpectralCube *cube = cube_from_values(1U, 2U, 4U, values);
    if (cube == NULL) {
        return EXIT_FAILURE;
    }
    const int status = spectral_cube_write_file(path, cube);
    spectral_cube_free(cube);
    return status == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

static int verify_manual_cli_fixture(const char *path) {
    SpectralCube *cube = spectral_cube_read_file(path);
    if (cube == NULL) {
        return EXIT_FAILURE;
    }
    int valid =
        cube->height == 1U && cube->width == 2U && cube->bands == 4U;
    size_t element_count = 0U;
    if (valid != 0 &&
        spectral_cube_element_count(cube, &element_count) == 0) {
        for (size_t index = 0U; index < element_count; ++index) {
            if (!isfinite(cube->data[index])) {
                valid = 0;
            }
        }
    } else {
        valid = 0;
    }
    spectral_cube_free(cube);
    return valid != 0 ? EXIT_SUCCESS : EXIT_FAILURE;
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

int main(int argc, char *argv[]) {
    if (argc == 3 && strcmp(argv[1], "--write-cli-fixture") == 0) {
        return write_manual_cli_fixture(argv[2]);
    }
    if (argc == 3 && strcmp(argv[1], "--verify-cli-fixture") == 0) {
        return verify_manual_cli_fixture(argv[2]);
    }
    if (argc != 1) {
        return EXIT_FAILURE;
    }

    test_allocation_and_copy_helpers();
    test_3d_helpers();
    test_circular_shifts();
    test_filtering_and_convolution();
    test_mirror_and_sampling();
    test_haar_dwt();
    test_multilevel_haar_dwt();
    test_direct_haar_thresholding();
    test_haar_detail_thresholding();
    test_haar_denoising();
    test_spectral_cube_model();
    test_spectral_cube_denoising();
    test_spectral_cube_file_io();
    test_spectral_denoise_cli();
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
