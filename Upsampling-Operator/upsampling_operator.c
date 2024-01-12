#include <stdio.h>
#include <stdlib.h>

#include "upsampling_operator.h"
#include "../Universal-Tools/universal_tools.h"

/* Created:       20.11.2023
   Last modified: 20.11.2023
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to perform umpsamling operation */
double ** upsampling_operator(double ** input, int rows, int columns, int scale) {

    int column_count_new = columns * scale;
    double ** upsampled = create_2d_array(rows, column_count_new);
    column_copying_worker(input, upsampled,  rows, columns, scale);
    
    return upsampled;
}
