#ifndef HIGH_PASS_UPSAMPLING_OPERATOR_H
#define HIGH_PASS_UPSAMPLING_OPERATOR_H

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double ** high_pass_upsampling_operator(double ** input, double * filter, int rows, int columns, int scale, int filter_length);

#endif  // HIGH_PASS_UPSAMPLING_OPERATOR_H
