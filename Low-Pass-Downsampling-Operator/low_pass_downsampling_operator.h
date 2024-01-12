#ifndef LOW_PASS_DOWNSAMPLING_OPERATOR_H
#define LOW_PASS_DOWNSAMPLING_OPERATOR_H

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double ** low_pass_downsampling_operator(double ** input, double * filter, int rows, int columns, int filter_length);

#endif  // LOW_PASS_DOWNSAMPLING_OPERATOR_H
