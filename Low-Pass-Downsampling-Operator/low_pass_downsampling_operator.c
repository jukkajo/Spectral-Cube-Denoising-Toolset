#include <stdio.h>
#include "../Low-Pass-Downsampling-Operator/low_pass_downsampling_operator.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Time-Reversed-Periodic-Convolution/time_reversed_periodic_convolution.h"

/* Created:       13.12.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* 2D low-pass downsampling operator */
double ** low_pass_downsampling_operator(double ** input, double * filter, int rows, int columns, int filter_length) {
                                     
    double ** result = time_reversed_periodic_convolution(input, filter, rows, columns, filter_length);
    
    /* Allocating memory and copying every second column */
    int new_column_count = columns / 2;
    double ** output = create_2d_array(rows, new_column_count);
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < new_column_count; j++) {
            output[i][j] = result[i][2 * j];
        }
    }
    
    return output;
}

