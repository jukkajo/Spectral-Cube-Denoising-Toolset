#include <stdio.h>
#include <stdlib.h>

#include "high_pass_upsampling_operator.h"
#include "../Mirror-Filter/mirror_filter.h"
#include "../Upsampling-Operator/upsampling_operator.h"
#include "../Circular-Right-Shift/circular_right_shift.h"
#include "../Periodic-Convolution/periodic_convolution.h"

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to perform high-pass upsampling */
double ** high_pass_upsampling_operator(double ** input, double * filter, int rows, int columns, int scale, int filter_length) {
  
    mirror_filter(filter, filter_length);
    
    double ** upsampled = upsampling_operator(input, rows, columns, scale);
    
    two_dim_circular_right_shift(upsampled, rows, (columns * scale));

    double ** result = periodic_convolution(filter, filter_length, upsampled, rows, (columns * scale));
    return result;

}
