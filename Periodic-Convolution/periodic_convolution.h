#ifndef TWO_DIM_CONVOLUTION_H
#define TWO_DIM_CONVOLUTION_H

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double ** two_dim_convolution(double * filter, int filter_length, double ** input, int rows, int columns);

#endif  // TWO_DIM_CONVOLUTION_H
