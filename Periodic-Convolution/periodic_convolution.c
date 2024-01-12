#include <stdio.h>
#include <stdlib.h>

#include "periodic_convolution.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Rational-Transfer-Function/rational_transfer_function.c"

/* Created:       20.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function for two-dimensional convolution, for
   two-scale transform. Filter is directly applied.
   Compatible with cases where:
   filter length < column count 
   and column count <= filter length
*/
double ** two_dim_convolution(double * filter, int filter_length, double ** input, int rows, int columns) {
    
    double ** padded_input;
    int new_column_count;
    
    /* If filter is smaller than column count, we pad the original input
       matrix with its own last columns.
       "periodic extension of the original data"
    */
    if (filter_length <= columns) {

        int start_column = columns - filter_length;
        int submatrix_columns = columns - start_column;
        
        /* Function to copy submatrix column-wise (all rows) from 2D array */
        double ** sub_matrix = column_wise_2d_array_submatrix_duplication(input, rows, columns, start_column, columns);
        
        /* Allocating memory for padded matrix */
        new_column_count = columns + submatrix_columns;
        
        /* Assigning padding (submatrix) to left side of the original input */
        padded_input = create_2d_array(rows, new_column_count);
        
        /* Copying first submatrix values */
        for (int i = 0; i < submatrix_columns; i++) {
            for (int j = 0; j < rows; j++) {
                padded_input[j][i] = sub_matrix[j][i];
            }
        }
        free_2d_array(sub_matrix, rows);

        /* Then, copying the original input matrix */
        for (int i = submatrix_columns; i < new_column_count; i++) {
            for (int j = 0; j < rows; j++) {
                padded_input[j][i] = input[j][i-submatrix_columns];   
            }
        }

    } else {
        /* Creating matrix of zeroes (padding), all rows and filter size amount of columns */
        double ** padding_array = create_2d_array(rows, filter_length);

        for (int i = 0; i < filter_length; i++) {
        
            int column_index = (filter_length * columns - filter_length + i) % columns;
            
            /* Column i (all rows) is initialized with column specified by column_index */
            for (int j = 0; j < rows; j++) {
                padding_array[j][i] = input[j][column_index];
            }

        }
        
        new_column_count = columns + filter_length - 1;
        padded_input = create_2d_array(rows, new_column_count);
        
        /* Copying first values from padding array */
        for (int i = 0; i < filter_length; i++) {
            for (int j = 0; j < rows; j++) {
                padded_input[j][i] = padding_array[j][i];
            }
        }
        free_2d_array(padding_array, rows);

        /* Then, copying the original input matrix.
           Padding end + 1 => starting index
        */
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                padded_input[j][i + filter_length] = input[j][i];
            }
        }

    }
    
    double ** convolved = rational_transfer_function(padded_input, rows, new_column_count, filter, filter_length);
    free_2d_array(padded_input, rows);
    
    /* Allocating memory for output array accordingly */
    double ** output = create_2d_array(rows, columns);
    
    /* Extracting all rows, output dimensions remains same as input's */
    for (int i = 0; i < rows; i++) {
        for (int j = filter_length; j < columns + filter_length; j++) {
            output[i][j - filter_length] = convolved[i][j];
        }
    }
    free_2d_array(convolved, rows);
    
    return output;
}
