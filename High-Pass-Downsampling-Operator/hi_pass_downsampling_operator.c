#include <stdio.h>
#include <stdlib.h>

#include "../High-Pass-Downsampling/high_pass_downsampling_operator.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Circular-Left-Shift/circular_left_shift.h"
#include "../Periodic-Convolution/periodic_convolution.h"
#include "../Mirror-Filter/mirror_filter.h"

/* Created:       20.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* 2D high-pass downsampling operator */
double ** hi_pass_downsampling_operator(double ** input, double * filter, int rows, int columns, int filter_length) {

    /* Implementing Nyquist theorem, function modifies values in filter's memory positions accordingly */
    mirror_filter(filter, filter_length);
    
    /* Aplying circular left shift on columns of matrix */
    two_dim_circular_left_shift(input, rows, columns);
    
    /* Performing (column-wise) convolution */
    double ** result = two_dim_convolution(filter, filter_length, input, rows, columns);

    /* Allocating memory for output */
    int new_column_count = columns / 2;
    double ** output = create_2d_array(rows, new_column_count);
    
    /* Copying every second column to output*/
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < new_column_count; j++) {
            output[i][j] = result[i][2 * j];
        }
    }
    
    /* Freeing allocated memory */
    free_2d_array(result, rows);
    
    return output;
}
