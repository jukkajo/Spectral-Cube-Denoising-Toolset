#ifndef MIRROR_FILTER_H
#define MIRROR_FILTER_H

#include <stddef.h>

/*
 * Return a caller-owned copy with alternating signs: output[k] is filter[k]
 * for even k and -filter[k] for odd k. The input is never modified. A null
 * filter or zero length returns NULL with errno set to EINVAL. Allocation
 * failure returns NULL; the successful result is caller-owned.
 */
double *mirror_filter(const double *filter, size_t filter_length);

#endif
