#ifndef DISCRETE_WAVELET_TRANSFORM_H
#define DISCRETE_WAVELET_TRANSFORM_H

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

double *** discrete_wavelet_transform(double *** input, double * filter, int filter_length, int rows, int columns, int pages, int scale);

#endif  // DISCRETE_WAVELET_TRANSFORM_H
