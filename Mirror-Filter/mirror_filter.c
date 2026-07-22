#include "mirror_filter.h"

#include <errno.h>

#include "../Universal-Tools/universal_tools.h"

double *mirror_filter(const double *filter, size_t filter_length) {
    if (filter == NULL || filter_length == 0U) {
        errno = EINVAL;
        return NULL;
    }

    double *mirrored = create_1d_array(filter_length);
    if (mirrored == NULL) {
        return NULL;
    }
    for (size_t index = 0U; index < filter_length; ++index) {
        mirrored[index] = (index % 2U == 0U) ? filter[index] : -filter[index];
    }
    return mirrored;
}
