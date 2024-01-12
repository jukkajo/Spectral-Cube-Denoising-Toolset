#include <stdio.h>
#include <stdlib.h>

#include "./Universal-Tools/universal_tools.h"
#include "./Discrete-Wavelet-Transform/discrete_wavelet_transform.h"
#include "./Wavelet-Coefficients/wavelet_coefficients.h"

/* main.c, for testing purposes, and to show example
   function call sequence
*/
int main(int argc, char *argv[]) {
    printf("T0");
    /* Empty 3D array of zeroes, imitates spectral cube. */
    
    double *** input = create_3d_zero_array(200, 200, 200);

    // TODO: DELETE
    double grow = 0.0;
    for (int i = 0; i < 200; i++) {
        for (int j = 0; j < 200; j++) {
            for (int k = 0; k < 200; k++) {

                input[j][k][i] = grow;
                grow += 0.0000003;
            }
        }

    }
    // TODO: DELETE

    printf("T1");
    int scale = 2;

    double *** transformed = discrete_wavelet_transform(input, coiflet_12, length_coiflet_12, 200, 200, 200, scale);
    
    //test
    double test[3][3][2] = {
        {{0.4344, 0.7954, 0.5675},
         {0.6021, 0.8763, 0.3268},
         {0.6385, 0.6908, 0.4215}},
    
        {{0.7421, 0.2763, 0.6803},
         {0.2739, 0.4034, 0.1942},
         {0.1032, 0.2200, 0.6108}}
    };
    int test_order[3] = {1,3,2};
    
    double *** tested = permute_3d(test, 3, 3, 2, test_order);
  
}
