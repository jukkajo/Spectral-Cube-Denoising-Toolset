#include <stdio.h>
#include <stdlib.h>

#include "time_reversed_periodic_convolution.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Rational-Transfer-Function/rational_transfer_function.c"

/* Created:       23.11.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to perform filtering by periodic convolution,
   it convolves matrix x with the time-reverse of filter f
   columnwise.
*/
double ** time_reversed_periodic_convolution(double ** input, double * filter, int rows, int columns, int filter_length) {

    double ** padded_input;
    int new_column_count;
    
    /* If filter is smaller than column count, we pad the original input
       matrix with its own last columns. A.k.a " circular or periodic extension
       of the original data"
    */
    if (filter_length <= columns) {

        int padding_start_column = 0;
        int padding_column_count = filter_length;
        /* Function to copy submatrix column-wise (all rows) from 2D array */
        double ** sub_matrix = column_wise_2d_array_submatrix_duplication(input, rows, columns, padding_start_column, padding_column_count);
        
        /* Allocating memory for padded matrix */
        new_column_count = columns + padding_column_count;
        padded_input = create_2d_array(rows, new_column_count);
        
        /* Copying first original input matrix */
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                padded_input[j][i] = input[j][i];
            }
        }
        
        /* Then, copying the submatrix values */

        for (int i = 0; i < (new_column_count - columns); i++) {
            for (int j = 0; j < rows; j++) {
               padded_input[j][columns + i] = sub_matrix[j][i];
            }
        }

        /* Now we have padded matrix */
    }  else {
        
        double ** padding_array = create_2d_array(rows, filter_length);
        
        for (int i = 0; i < filter_length; i++) {
           int column_index = i % columns;
           
           /* Column i (all rows) is initialized with column specified by column_index */
           for (int j = 0; j < rows; j++) {
               padding_array[j][i] = input[j][column_index];
           }
        }
        
        new_column_count = columns + filter_length;
        padded_input = create_2d_array(rows, new_column_count);
        
        /* Copying first original input matrix */
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < rows; j++) {
                padded_input[j][i] = input[j][i];
            }
        }

        /* Then, copying the submatrix values */

        for (int i = 0; i < (new_column_count - columns); i++) {
            for (int j = 0; j < rows; j++) {
               padded_input[j][columns + i] = padding_array[j][i];
            }
        }

    }
    
    double * reversed_filter = reverse_1d_array(filter, filter_length);
    double ** convolved = rational_transfer_function(padded_input, rows, new_column_count, reversed_filter, filter_length);
    double ** output = create_2d_array(rows, columns);

    /* Extracting all rows, in defined column range of:
       (filter_length - 1) to new_column_count
       Therefore, output dimensions are same as input's
    */

    for (int i = 0; i < rows; i++) {
        for (int j = filter_length; j < new_column_count; j++) {
            output[i][j - filter_length] = convolved[i][j-1];
        }
    }
    
    return output;
}
