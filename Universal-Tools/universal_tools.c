#include <stdio.h>
#include <stdlib.h>

#include "./universal_tools.h"

/* Created:       20.11.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to free dynamically allocated memory, for 3D structures  */
void free_3d_array(double ***array, int rows, int columns) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}

/* Function to free dynamically allocated memory, for 2D structures */
void free_2d_array(double ** array, int rows) {
    for (int i = 0; i < rows; i++) {
        free(array[i]);
    }
    free(array);
}

/* Function to create and initialize a 1D array as zeros */
double * create_1d_array(int length) {
    double * array;
    array = (double *)malloc(length * sizeof(double));
    return array;
}

/* Function to create and initialize a 2D array as zeros */
double ** create_2d_array(int rows, int columns) {
    double ** array = (double**)malloc(rows * sizeof(double *));

    for (int i = 0; i < rows; i++) {
        array[i] = (double *)malloc(columns * sizeof(double));

        for (int j = 0; j < columns; j++) {
            array[i][j] = 0.0;
        }
    }
    return array;
}

/* Function to create and initialize a 3D array as zeros */
double *** create_3d_zero_array(int rows, int columns, int pages) {
    double *** three_dim_array = (double ***)malloc(rows * sizeof(double **));

    for (int i = 0; i < rows; ++i) {
        three_dim_array[i] = (double **)malloc(columns * sizeof(double *));
        
        for (int j = 0; j < columns; ++j) {
            three_dim_array[i][j] = (double *)malloc(pages * sizeof(double));

            for (int k = 0; k < pages; ++k) {
                three_dim_array[i][j][k] = 0.0;
            }
        }
    }

    return three_dim_array;
}

/* Function to flatten 3D array to 2D */
double ** flatten_3d_to_2d(double *** input, int rows, int columns, int pages) {

    double ** flattened = create_2d_array((rows * columns), pages);
    for (int i = 0; i < pages; i++) {
        int counter = 0;
        int a = 0;
        int b = 0;
        
        for (int j = 0; j < columns; j++) {
            for (int k = 0; k < rows; k++) {
                flattened[counter][i] = input[k][j][i]; 
                counter++;
                b = k;
            }
            a = j;
        }
    }
    return flattened;
}

/* Function to copy submatrix column-wise (all rows) from 2D array */
double ** column_wise_2d_array_submatrix_duplication(double ** source, int rows_original, int columns_original, int column_range_start, int column_range_end) {

    int submatrix_columns = column_range_end - column_range_start;
    double ** sub_matrix = (double **)malloc(rows_original * sizeof(double *));
    for (int i = 0; i < rows_original; i++) {
        sub_matrix[i] = (double *)malloc(submatrix_columns * sizeof(double));
        for (int j = 0; j < submatrix_columns; j++) {
            sub_matrix[i][j] = source[i][column_range_start + j];
        }
    }

    return sub_matrix;
}

/* Function to copy desired 2D submatrices within given range (depth-wise) */
double *** depth_wise_3d_submatrix_duplication(double *** input, int rows, int columns, int depth_range[2]) {
 
    int pages = depth_range[1] - depth_range[0] + 1;
    
    /* Allocating memory for output */
    double *** submatrix = create_3d_zero_array(rows, columns, pages);
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            int iterator = 0;
            for (int k = depth_range[0]; k < depth_range[1]; k++) {
                submatrix[i][j][iterator] = input[i][j][k];
                iterator++;
            }
        }
    }
    
    return submatrix;
}

/* Function to copy submatrix column-wise (all rows) from 2D array
   Build to be universal tool to some extent
*/
void column_copying_worker(double **input, double **output, int rows, int columns, int step) {
    for (int j = 0; j < columns; j++) {
        for (int i = 0; i < rows; i++) {
            output[i][(j * step)] = input[i][j];
        }
    }
}

/* Function to reverse a 1D array */
double * reverse_1d_array(double * array, int length) {
    double * reversed = (double *)malloc(length * sizeof(double));

    int start = 0;
    int end = length - 1;

    while (start < length) {
        reversed[start] = array[end];
        start++;
        end--;
    }

    return reversed;
}

/* Function to reshape a matrix:
   from 2d to 3d or from 2d to 3d
   
   returns structs containing either dimensioned matrix and null-pointer
*/
struct Reshaped reshape(double ** input_1, double *** input_2, int bottom_z[2], int top_z[2], int dimension_variable, int rows, int columns, int pages) {
    struct Reshaped reshaped;
    
    int iter1, iter2, iter3, counter;
    iter1 = 0; iter2 = 0; iter3 = 0; counter = 0;
    
    if (input_1 != NULL && input_2 == NULL) {
    
        reshaped.matrix_3d = create_3d_zero_array(rows, columns, pages);
        printf("Dim-var: %d\n", dimension_variable);
        while (iter3 < dimension_variable) {
            while (iter1 < columns) {
               //printf("Iter2: %d ---- Iter1: %d ---- Counter: %d\n", iter2, iter1, counter);
               reshaped.matrix_3d[iter2][iter1][counter] = input_1[iter3][iter1];
               iter1++;
            }
            /* ---------------- */
            //printf("Counter: %d\n", iter2);
            if (iter2 < rows-1) {
                iter2++;
            } else {
                iter2 = 0;
                 
                counter++;
            }
            iter3++;
            iter1 = 0;
            /* ---------------- */
        }
        
    } else if (input_2 != NULL && input_1 == NULL) {
        reshaped.matrix_2d = create_2d_array(rows, columns);
        /* TODO: Implement if needed, remove otherwise */
    }

    return reshaped;
}

/* Function to rearrange dimensions of 3D matrix */
double *** permute_3d(double *** input, int rows,  int columns, int pages, int new_order[3]) {
    double *** permuted;
    
    /* TODO: add more usecases if needed */
    switch (new_order[0]) {
        case 1:
            /* 1,3,2 */
            permuted = create_3d_zero_array(rows, pages, columns);
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < columns; j++) {
                    for (int k = 0; k < pages; k++) {
                        if(k < ) {
                        } else {
                        
                        }
                        permuted[i][j][k] = input[i][j][k];
                        
                    }
                }
            }
            break;

        case 2:
            /* 2,1,3 */
            permuted = create_3d_zero_array(columns, rows, pages);
            // Add code to rearrange dimensions for case 2
            break;

        // Add more cases if needed

    }
    
    return permuted;
}


