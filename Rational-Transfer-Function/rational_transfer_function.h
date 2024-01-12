#ifndef RATIONAL_TRANSFER_FUNCTION_H
#define RATIONAL_TRANSFER_FUNCTION_H

/* Created:       23.11.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double ** rational_transfer_function(double ** input, int rows, int columns, double * filter, int filter_length);

#endif  // RATIONAL_TRANSFER_FUNCTION_H
