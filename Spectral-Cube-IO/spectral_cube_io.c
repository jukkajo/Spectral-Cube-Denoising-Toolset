#define _POSIX_C_SOURCE 200809L

#include "spectral_cube_io.h"

#include <errno.h>
#include <float.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TEMPORARY_FILE_ATTEMPTS 100U
#define TEMPORARY_SUFFIX_CAPACITY 32U

static int validate_binary64_host(void) {
    if (CHAR_BIT != 8 || sizeof(double) != 8U || FLT_RADIX != 2 ||
        DBL_MANT_DIG != 53 || DBL_MAX_EXP != 1024 ||
        DBL_MIN_EXP != -1021) {
        errno = ENOTSUP;
        return -1;
    }

    const double one = 1.0;
    uint64_t one_bits = 0U;
    memcpy(&one_bits, &one, sizeof(one_bits));
    if (one_bits != UINT64_C(0x3ff0000000000000)) {
        errno = ENOTSUP;
        return -1;
    }
    return 0;
}

static void encode_uint64_le(uint64_t value, unsigned char bytes[8]) {
    for (size_t index = 0U; index < 8U; ++index) {
        bytes[index] = (unsigned char)(value & UINT64_C(0xff));
        value >>= 8U;
    }
}

static uint64_t decode_uint64_le(const unsigned char bytes[8]) {
    uint64_t value = 0U;
    for (size_t index = 0U; index < 8U; ++index) {
        const unsigned int shift = (unsigned int)(8U * index);
        value |= (uint64_t)bytes[index] << shift;
    }
    return value;
}

static void encode_double_le(double value, unsigned char bytes[8]) {
    uint64_t bits = 0U;
    memcpy(&bits, &value, sizeof(bits));
    encode_uint64_le(bits, bytes);
}

static double decode_double_le(const unsigned char bytes[8]) {
    const uint64_t bits = decode_uint64_le(bytes);
    double value = 0.0;
    memcpy(&value, &bits, sizeof(value));
    return value;
}

static void set_stream_error(void) {
    if (errno == 0) {
        errno = EIO;
    }
}

static int read_exact(FILE *stream, unsigned char *buffer, size_t length) {
    size_t offset = 0U;
    while (offset < length) {
        errno = 0;
        const size_t amount = fread(buffer + offset, 1U, length - offset, stream);
        if (amount == 0U) {
            if (ferror(stream) != 0) {
                set_stream_error();
            } else {
                errno = EINVAL;
            }
            return -1;
        }
        offset += amount;
    }
    return 0;
}

static int write_exact(
    FILE *stream,
    const unsigned char *buffer,
    size_t length
) {
    size_t offset = 0U;
    while (offset < length) {
        errno = 0;
        const size_t amount =
            fwrite(buffer + offset, 1U, length - offset, stream);
        if (amount == 0U) {
            set_stream_error();
            return -1;
        }
        offset += amount;
    }
    return 0;
}

static int uint64_to_size(uint64_t value, size_t *converted) {
    const size_t result = (size_t)value;
    if ((uint64_t)result != value) {
        errno = EOVERFLOW;
        return -1;
    }
    *converted = result;
    return 0;
}

static int size_to_uint64(size_t value, uint64_t *converted) {
    const uint64_t result = (uint64_t)value;
    if ((size_t)result != value) {
        errno = EOVERFLOW;
        return -1;
    }
    *converted = result;
    return 0;
}

static SpectralCube *fail_read(FILE *stream, SpectralCube *cube, int error) {
    if (stream != NULL) {
        (void)fclose(stream);
    }
    spectral_cube_free(cube);
    errno = error;
    return NULL;
}

SpectralCube *spectral_cube_read_file(const char *path) {
    if (path == NULL || path[0] == '\0') {
        errno = EINVAL;
        return NULL;
    }
    if (validate_binary64_host() != 0) {
        return NULL;
    }

    FILE *stream = fopen(path, "rb");
    if (stream == NULL) {
        return NULL;
    }

    unsigned char header[SPECTRAL_CUBE_FILE_HEADER_SIZE];
    if (read_exact(stream, header, sizeof(header)) != 0) {
        const int saved_errno = errno;
        return fail_read(stream, NULL, saved_errno);
    }
    if (memcmp(header, SPECTRAL_CUBE_FILE_MAGIC, 8U) != 0) {
        return fail_read(stream, NULL, EINVAL);
    }

    size_t height = 0U;
    size_t width = 0U;
    size_t bands = 0U;
    if (uint64_to_size(decode_uint64_le(header + 8U), &height) != 0 ||
        uint64_to_size(decode_uint64_le(header + 16U), &width) != 0 ||
        uint64_to_size(decode_uint64_le(header + 24U), &bands) != 0) {
        const int saved_errno = errno;
        return fail_read(stream, NULL, saved_errno);
    }

    SpectralCube *cube = spectral_cube_create(height, width, bands);
    if (cube == NULL) {
        const int saved_errno = errno;
        return fail_read(stream, NULL, saved_errno);
    }

    size_t element_count = 0U;
    if (spectral_cube_element_count(cube, &element_count) != 0) {
        const int saved_errno = errno;
        return fail_read(stream, cube, saved_errno);
    }

    unsigned char encoded_value[8];
    for (size_t index = 0U; index < element_count; ++index) {
        if (read_exact(stream, encoded_value, sizeof(encoded_value)) != 0) {
            const int saved_errno = errno;
            return fail_read(stream, cube, saved_errno);
        }
        cube->data[index] = decode_double_le(encoded_value);
    }

    errno = 0;
    const int trailing = fgetc(stream);
    if (trailing != EOF) {
        return fail_read(stream, cube, EINVAL);
    }
    if (ferror(stream) != 0) {
        set_stream_error();
        const int saved_errno = errno;
        return fail_read(stream, cube, saved_errno);
    }
    if (fclose(stream) != 0) {
        const int saved_errno = errno == 0 ? EIO : errno;
        spectral_cube_free(cube);
        errno = saved_errno;
        return NULL;
    }

    return cube;
}

static FILE *open_temporary_file(
    const char *destination_path,
    char **temporary_path
) {
    const size_t destination_length = strlen(destination_path);
    if (destination_length > SIZE_MAX - TEMPORARY_SUFFIX_CAPACITY) {
        errno = EOVERFLOW;
        return NULL;
    }

    const size_t capacity = destination_length + TEMPORARY_SUFFIX_CAPACITY;
    char *path = malloc(capacity);
    if (path == NULL) {
        errno = ENOMEM;
        return NULL;
    }

    for (unsigned int attempt = 0U;
         attempt < TEMPORARY_FILE_ATTEMPTS;
         ++attempt) {
        const int written = snprintf(
            path,
            capacity,
            "%s.tmp.%u",
            destination_path,
            attempt
        );
        if (written < 0 || (size_t)written >= capacity) {
            free(path);
            errno = EOVERFLOW;
            return NULL;
        }

        FILE *stream = fopen(path, "wbx");
        if (stream != NULL) {
            *temporary_path = path;
            return stream;
        }
        if (errno != EEXIST) {
            const int saved_errno = errno;
            free(path);
            errno = saved_errno;
            return NULL;
        }
    }

    free(path);
    errno = EEXIST;
    return NULL;
}

static int fail_write(
    FILE *stream,
    char *temporary_path,
    int error
) {
    if (stream != NULL) {
        (void)fclose(stream);
    }
    if (temporary_path != NULL) {
        (void)unlink(temporary_path);
    }
    free(temporary_path);
    errno = error;
    return -1;
}

static int publish_temporary_file(
    const char *temporary_path,
    const char *destination_path
) {
    if (link(temporary_path, destination_path) != 0) {
        return -1;
    }

    errno = 0;
    if (unlink(temporary_path) != 0) {
        const int saved_errno = errno == 0 ? EIO : errno;
        (void)unlink(destination_path);
        errno = saved_errno;
        return -1;
    }
    return 0;
}

int spectral_cube_write_file(const char *path, const SpectralCube *cube) {
    if (path == NULL || path[0] == '\0') {
        errno = EINVAL;
        return -1;
    }
    if (validate_binary64_host() != 0) {
        return -1;
    }

    size_t element_count = 0U;
    if (spectral_cube_element_count(cube, &element_count) != 0) {
        return -1;
    }

    uint64_t height = 0U;
    uint64_t width = 0U;
    uint64_t bands = 0U;
    if (size_to_uint64(cube->height, &height) != 0 ||
        size_to_uint64(cube->width, &width) != 0 ||
        size_to_uint64(cube->bands, &bands) != 0) {
        return -1;
    }

    char *temporary_path = NULL;
    FILE *stream = open_temporary_file(path, &temporary_path);
    if (stream == NULL) {
        return -1;
    }

    unsigned char header[SPECTRAL_CUBE_FILE_HEADER_SIZE] = {0U};
    memcpy(header, SPECTRAL_CUBE_FILE_MAGIC, 8U);
    encode_uint64_le(height, header + 8U);
    encode_uint64_le(width, header + 16U);
    encode_uint64_le(bands, header + 24U);
    if (write_exact(stream, header, sizeof(header)) != 0) {
        const int saved_errno = errno;
        return fail_write(stream, temporary_path, saved_errno);
    }

    unsigned char encoded_value[8];
    for (size_t index = 0U; index < element_count; ++index) {
        encode_double_le(cube->data[index], encoded_value);
        if (write_exact(stream, encoded_value, sizeof(encoded_value)) != 0) {
            const int saved_errno = errno;
            return fail_write(stream, temporary_path, saved_errno);
        }
    }

    errno = 0;
    if (fflush(stream) != 0) {
        const int saved_errno = errno == 0 ? EIO : errno;
        return fail_write(stream, temporary_path, saved_errno);
    }
    errno = 0;
    if (fclose(stream) != 0) {
        const int saved_errno = errno == 0 ? EIO : errno;
        return fail_write(NULL, temporary_path, saved_errno);
    }
    stream = NULL;

    if (publish_temporary_file(temporary_path, path) != 0) {
        const int saved_errno = errno;
        return fail_write(NULL, temporary_path, saved_errno);
    }

    free(temporary_path);
    return 0;
}
