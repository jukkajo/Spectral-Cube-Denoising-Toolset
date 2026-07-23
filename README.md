# Spectral Cube Denoising Toolset

A C11 library and command-line program for applying multilevel Haar
thresholding independently to each pixel's spectral vector.

## Supported v0.1.0 scope

- One-level and multilevel normalized Haar transforms and reconstruction
- Even- and odd-length one-dimensional signals
- Hard and soft thresholding with user-provided thresholds
- Contiguous pixel-major `SpectralCube` storage
- Independent per-pixel spectral-axis denoising
- Version-1 `SCUBE001` file input and output
- The `spectral-denoise` file-to-file command
- Atomic POSIX no-overwrite output publication

## Requirements

- A C11-compatible compiler
- GNU Make
- A POSIX system with hard-link support
- The standard C and math libraries
- An 8-bit-byte, IEEE-754 binary64-compatible `double` representation for
  `SCUBE001` input and output

## Build and test

Build `libspectral_operators.a` and `spectral-denoise`:

```bash
make
```

Run the deterministic test suite:

```bash
make test
```

Run the suite with AddressSanitizer and UndefinedBehaviorSanitizer:

```bash
make sanitize
```

Remove generated build files:

```bash
make clean
```

## Command-line use

All processing options are required:

```bash
./spectral-denoise \
    --input noisy.scube \
    --output denoised.scube \
    --levels 2 \
    --threshold 0.75 \
    --mode soft
```

`--help` and `--version` are also available. Levels must be valid for the
input band count. Thresholds must be finite and nonnegative.

Hard thresholding uses:

```text
abs(value) < threshold  -> 0
otherwise               -> value
```

Values exactly equal to the threshold remain unchanged. Soft thresholding
uses:

```text
sign(value) * max(abs(value) - threshold, 0)
```

Values exactly equal to the threshold become zero. Approximation coefficients
are never thresholded.

## SCUBE001 format

```text
Offset  Size  Field
0       8     ASCII magic bytes SCUBE001
8       8     height as unsigned 64-bit little-endian
16      8     width as unsigned 64-bit little-endian
24      8     bands as unsigned 64-bit little-endian
32      ...   IEEE-754 binary64 values in little-endian order
```

Payload values use pixel-major order:

```text
((row * width) + column) * bands + band
```

Readers require the exact header and payload length. Writers complete and
close a temporary sibling before atomically publishing it with `link()`.
An existing destination, including a symbolic link, fails with `EEXIST` and
is never replaced. The temporary name is removed after success or failure.

## Current limitations

- Haar is the only implemented wavelet family.
- Thresholds are supplied by the user; there is no automatic noise or
  threshold estimation.
- Processing is limited to the spectral axis. Spatial and full 3D transforms
  are not implemented.
- `SCUBE001` carries dimensions and values only; it has no metadata.
- ENVI, FITS, HDF5, streaming, parallel processing, and GPU processing are not
  supported.
- The deterministic MSE regression validates only its recorded test fixture;
  it is not a general scientific performance claim.

## License status

No standalone license has been selected for this repository. The repository
owner must choose and add a license before creating a public release tag.
