#include <stdio.h>
#include <stdlib.h>

#include "mirror_filter.h"

/* Created:       23.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to shift 1D signal containg DC frequency content, to Nyquist frequency */
void mirror_filter(double * filter, int filter_length) {
    for (int i = 0; i < filter_length; i++) {
        filter[i] *= (i % 2 == 1) ? -1.0 : 1.0;
    }
}
