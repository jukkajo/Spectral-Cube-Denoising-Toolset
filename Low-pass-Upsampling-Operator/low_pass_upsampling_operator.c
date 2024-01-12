#include <stdio.h>
#include <stdlib.h>

#include "../Low-Pass-Upsampling-Operator/low_pass_upsampling_operator.h"
#include "../Universal-Tools/universal_tools.h"
#include "../Periodic-Convolution/periodic_convolution.h"

/* Created:       20.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to perform low-pass upsampling */
double ** low_pass_upsampling_operator(double ** input, double * filter, int rows, int columns, int scale, int filter_length) {

    double ** upsampled = upsampling_operator(input, rows, columns, scale);
    double ** result = two_dim_convolution(filter, filter_length, upsampled, rows, (columns * scale));
    free_2d_array(upsampled, rows);

    return result;

}
