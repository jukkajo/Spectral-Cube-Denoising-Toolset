#ifndef HIGH_PASS_DOWNSAMPLING_OPERATOR_H
#define HIGH_PASS_DOWNSAMPLING_OPERATOR_H

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double ** hi_pass_downsampling_operator(double ** input, double * filter, int rows, int columns, int filter_length);

#endif  // HIGH_PASS_DOWNSAMPLING_OPERATOR_H
