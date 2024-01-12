#include <stdio.h>
#include <stdlib.h>

#ifndef UNIVERSAL_TOOLS_H
#define UNIVERSAL_TOOLS_H

/* Created:       20.11.2023
   Last modified: 12.01.2024
   @ Jukka J jajoutzs@jyu.fi
*/

struct Reshaped {
    double ** matrix_2d;
    double *** matrix_3d;
};

void free_3d_array(double ***array, int rows, int columns);

void free_2d_array(double ** array, int rows);

double * create_1d_array(int length);

double ** create_2d_array(int rows, int columns);

double *** create_3d_zero_array(int rows, int columns, int pages);

double ** flatten_3d_to_2d(double *** input, int rows, int columns, int pages);

double ** column_wise_2d_array_submatrix_duplication(double ** source, int rows_original, int columns_original, int column_range_start, int column_range_end);

double *** depth_wise_3d_submatrix_duplication(double *** input, int rows, int columns, int depth_range[2]);
 
void column_copying_worker(double **input, double **output, int rows, int columns, int step);

double * reverse_1d_array(double * array, int length);

struct Reshaped reshape(double ** input_1, double *** input_2, int bottom_z[2], int top_z[2], int dimension_variable, int rows, int columns, int pages);

double *** permute_3d(double *** input, int rows,  int columns, int pages, int new_order[3]);

#endif // UNIVERSAL_TOOLS_H

