#ifndef SPECTRAL_CUBE_IO_H
#define SPECTRAL_CUBE_IO_H

#include "../Spectral-Cube/spectral_cube.h"

#define SPECTRAL_CUBE_FILE_MAGIC "SCUBE001"
#define SPECTRAL_CUBE_FILE_HEADER_SIZE 32U

/*
 * SCUBE001 version-1 binary layout
 * --------------------------------
 * Offset  Size  Field
 *      0     8  ASCII magic bytes "SCUBE001"
 *      8     8  height as an unsigned 64-bit little-endian integer
 *     16     8  width as an unsigned 64-bit little-endian integer
 *     24     8  bands as an unsigned 64-bit little-endian integer
 *     32   ...  IEEE-754 binary64 values as little-endian bytes
 *
 * Payload values follow SpectralCube's pixel-major order:
 *
 *   ((row * width) + column) * bands + band
 *
 * Only this version is supported. Native C structs are never written. File
 * operations first verify that the host double is an 8-byte IEEE-754 binary64
 * representation.
 */

/*
 * Read one complete SCUBE001 file into a caller-owned SpectralCube. Invalid
 * magic, metadata, truncation, or trailing bytes return NULL with errno set
 * to EINVAL. Dimension conversion or allocation arithmetic overflow uses
 * EOVERFLOW. Host incompatibility uses ENOTSUP. File-system errors are
 * preserved; stream errors without a more specific errno use EIO. No partial
 * cube is returned.
 */
SpectralCube *spectral_cube_read_file(const char *path);

/*
 * Write cube to a temporary sibling file, flush and close it, then publish it
 * without replacement using POSIX link() followed by unlink() of the
 * temporary name. Existing destination directory entries, including symbolic
 * links, are rejected atomically with EEXIST. Invalid arguments or cube
 * metadata use EINVAL, arithmetic conversion uses EOVERFLOW, and host
 * incompatibility uses ENOTSUP. File-system errors are preserved with EIO as
 * the stream-error fallback. Any temporary file created by this call is
 * removed after failure, and source is never modified.
 */
int spectral_cube_write_file(const char *path, const SpectralCube *cube);

#endif
