#ifndef TIME_REVERSED_PERIODIC_CONVOLUTION_H
#define TIME_REVERSED_PERIODIC_CONVOLUTION_H

/* Created:       23.11.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double ** time_reversed_periodic_convolution(double ** input, double * filter, int rows, int columns, int filter_length);

#endif  // TIME_REVERSED_PERIODIC_CONVOLUTION_H
