#ifndef CIRCULAR_RIGHT_SHIFT_H
#define CIRCULAR_RIGHT_SHIFT_H

#include <stddef.h>

/*
 * Rotate all columns right in place. Rows and columns must be positive and all
 * row pointers must be non-null. Returns 0 on success, or -1 with errno set.
 * A one-column matrix is a successful no-op.
 */
int circular_right_shift(double **input, size_t rows, size_t columns);

#endif
