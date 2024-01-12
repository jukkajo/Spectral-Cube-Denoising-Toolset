#include <stdio.h>
#include <stdlib.h>

#include "discrete_wavelet_transform.h"
#include "../Universal-Tools/universal_tools.h"
#include "../High-Pass-Downsampling-Operator/hi_pass_downsampling_operator.h"
#include "../Low-Pass-Downsampling-Operator/low_pass_downsampling_operator.h"

/* Created:       27.11.2023
   Last modified: 11.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

/* Function to compute wavelet transform for 3D-array */
double *** discrete_wavelet_transform(double *** input, double * filter, int filter_length, int rows, int columns, int pages, int scale) {
    //TODO:delete line below
    scale = 1; 
    int dimension_variable;
    
    printf("T2");
    for (int i = 0; i < scale; i++) {
        /* Ranges: */
        int top_x[2] = {(rows / 2), (rows - 1)};
        int top_y[2] = {(columns / 2), (columns - 1)};
        int top_z[2] = {(pages / 2), (pages - 1)};
        /*-------------------------------------*/
        int bottom_x[2] = {0, (rows / 2 - 1)};
        int bottom_y[2] = {0, (columns / 2 - 1)};
        int bottom_z[2] = {0, (pages / 2 - 1)};
        /*-------------------------------------*/
        printf("T3");
        int depth_range[2] = {0, pages};
        /* Copying submatrix. With the first iteration, it will copy entire input */
        double *** sub_matrix = depth_wise_3d_submatrix_duplication(input, rows, columns, depth_range);
        printf("T4");
        /* Creating 2D matrix */
        dimension_variable = rows * columns;
        printf("T5");

        // TODO: DELETE
        /*
        for (int i = 0; i < pages; i++) {
            printf("Page %d:\n", i + 1);

            for (int j = 0; j < rows; j++) {
                for (int k = 0; k < columns; k++) {
                    printf("%.7f\t", input[j][k][i]);
                }

                printf("\n");
            }

            printf("\n");
        }
        */
        // TODO: DELETE
        
        
        dimension_variable = rows * columns;
        double ** flattened = flatten_3d_to_2d(sub_matrix, rows, columns, pages);
        printf("T6");
        double ** hi_pass_downsampling_operator(flattened, filter, rows, columns, filter_length); 
        printf("T7");
        double ** low_pass_downsampling_operator(flattened, filter, rows, columns, filter_length);
        printf("T8");
        /*
        nx   ny     nz
        rows column pages
                  \ /
        ixy:       *
        
        */
        
        /* Within range ixy, take slices:
           
           bottom_z -> ... -> top_z
           
           reshape to rows column pages
        */
        printf("T9");
        /* ----------------- */
        struct Reshaped temporary = reshape(flattened, NULL, bottom_z, top_z, dimension_variable, rows, columns, pages);
        double *** reshaped = temporary.matrix_3d;

        // TODO: DELETE
        /*
        for (int i = 0; i < pages; i++) {
            printf("Page %d:\n", i + 1);

            for (int j = 0; j < rows; j++) {
                for (int k = 0; k < columns; k++) {
                    printf("%.7f\t", reshaped[j][k][i]);
                }

                printf("\n");
            }

            printf("\n");
        }
        */
        // TODO: DELETE
        
        int copy_all_pages[2] = {0, columns};
        double *** copy_of_reshaped = depth_wise_3d_submatrix_duplication(reshaped, rows, columns, copy_all_pages);
        
        
        /*-------------------*/
        rows = rows / 2;
        columns = columns / 2;
        pages = pages / 2;
    }
    
    
}

