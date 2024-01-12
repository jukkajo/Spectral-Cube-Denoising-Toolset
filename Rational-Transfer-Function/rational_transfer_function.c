#include <stdio.h>
#include <stdlib.h>

#include "rational_transfer_function.h"
#include "../Universal-Tools/universal_tools.h"

/* Created:       23.11.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to filter the input data, by using method of rational transfer function */
double ** rational_transfer_function(double ** input, int rows, int columns, double * filter, int filter_length) {

    /* Matrix to store the convolution result */
    double ** output = create_2d_array(rows, columns);

    /* Performing convolution */
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < columns; j++) {
            /* Assigning to zero, so we can sum to it */
            output[i][j] = 0.0;

            for (int k = 0; k < filter_length; k++) {
                /* Ensuring that we don't go out of bounds */
                if (j - k >= 0) {
                    output[i][j] += input[i][j - k] * filter[k];
                }
            }
        }
    }

    return output;

}
