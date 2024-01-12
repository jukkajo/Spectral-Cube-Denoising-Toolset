#include <stdio.h>
#include <stdlib.h>

#include "circular_left_shift.h"
#include "../Universal-Tools/universal_tools.h"

/* Created:       20.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to perform circular left shift on the columns of 2D input-signal */
void circular_left_shift(double ** input, int rows, int columns) {
    
    /* Temporary array to hold the values of the first column */
    double * temp = create_1d_array(rows);

    /* Storing the values of the first column */
    for (int i = 0; i < rows; i++) {
        temp[i] = input[i][0];
    }

    /* Shifting columns to left */
    for (int j = 0; j < columns - 1; j++) {
        for (int i = 0; i < rows; i++) {
            input[i][j] = input[i][j + 1];
        }
    }

    /* Moving values from the temporary array to the rightmost column */
    for (int i = 0; i < rows; i++) {
        input[i][columns - 1] = temp[i];
    }
    free(temp);
}
