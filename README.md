# Spectral Cube Denoising Toolset

A C11 numerical operator library and focused command-line workflow for Haar denoising along the spectral axis of a contiguous spectral cube.

The current implementation provides tested signal-processing primitives, multidimensional array utilities, an immutable wavelet coefficient catalogue, reversible one-level and multilevel one-dimensional Haar transforms, user-directed coefficient thresholding, independent spectral-axis denoising, the versioned `SCUBE001` binary cube format, and the `spectral-denoise` executable. Non-Haar transforms, automatic threshold selection, metadata, and other scientific file formats are not implemented.

## Status

The numerical operator layer is implemented and verified.

- The project builds without compiler warnings.
- All operator implementations link exactly once into a static library.
- The test suite contains 5,864 deterministic checks.
- AddressSanitizer and UndefinedBehaviorSanitizer checks pass.
- One-level and multilevel normalized Haar forward and inverse transforms support even and odd one-dimensional signals.
- Hard and soft thresholding support user-selected Haar detail-level ranges.
- Contiguous pixel-major spectral cubes support independent per-pixel spectral denoising.
- `SCUBE001` files round-trip dimensions and binary64 values exactly.
- `spectral-denoise` provides one file-to-file Haar spectral-axis workflow.
- The reserved generic 3D transform API remains unsupported and returns `ENOTSUP` for valid requests.

## Build

Requirements:

- A C11-compatible compiler
- GNU Make
- Standard C library
- Math library

Build the library and executable:

```bash
make
```

This produces:

```text
libspectral_operators.a
spectral-denoise
```

Remove generated build files:

```bash
make clean
```

## Test

Run the operator test suite:

```bash
make test
```

Run the tests with AddressSanitizer and UndefinedBehaviorSanitizer:

```bash
make sanitize
```

LeakSanitizer can be enabled on compatible systems:

```bash
make sanitize LEAK_CHECK=1
```

Some restricted or traced environments do not support LeakSanitizer.

## Spectral-cube file format

Only `SCUBE001` version 1 is supported. Integers and floating-point values are encoded explicitly in little-endian byte order; native C structs are not written.

```text
Offset  Size  Field
0       8     ASCII magic bytes SCUBE001
8       8     height as unsigned 64-bit little-endian
16      8     width as unsigned 64-bit little-endian
24      8     bands as unsigned 64-bit little-endian
32      ...   IEEE-754 binary64 values in little-endian order
```

Payload values use the same pixel-major order as `SpectralCube`:

```text
((row * width) + column) * bands + band
```

Readers require an exact header and payload length and reject trailing bytes. Writers publish a completed temporary sibling with an atomic POSIX hard link, so an existing destination is never replaced; the temporary name is removed after success or failure.

## Command-line workflow

The executable applies the existing multilevel Haar denoising operation independently to each pixel's spectral vector:

```bash
./spectral-denoise \
    --input noisy.scube \
    --output denoised.scube \
    --levels 2 \
    --threshold 0.75 \
    --mode soft
```

All five processing options are required. `--mode` accepts only `hard` or `soft`; `--help` and `--version` are also available.

## Implemented components

| Component | Purpose |
|---|---|
| Circular left and right shift | Circularly reorders one-dimensional signals |
| Periodic convolution | Applies convolution with periodic boundary indexing |
| Time-reversed periodic convolution | Applies the documented reversed-filter convolution convention |
| Rational transfer function | Applies causal finite impulse response filtering |
| Upsampling | Inserts zero-valued samples between input values |
| Low-pass downsampling | Filters and retains samples at even indices |
| High-pass downsampling | Applies a mirrored high-pass filter before downsampling |
| Low-pass upsampling | Upsamples and applies low-pass filtering |
| High-pass upsampling | Upsamples and applies high-pass filtering |
| Mirror filter | Creates a quadrature-style mirrored filter without modifying the input |
| Universal tools | Allocation, copying, slicing, reversal, flattening, reshaping, and permutation utilities |
| Wavelet coefficients | Immutable coefficient sets for supported wavelet families |
| Discrete wavelet transform | Reversible one-level and multilevel 1D Haar transforms with per-level odd-length metadata |
| Haar denoising | User-directed hard or soft detail thresholding followed by reconstruction |
| Spectral cube | Contiguous pixel-major storage and independent spectral-axis Haar denoising |
| Spectral cube I/O | Exact version-1 `SCUBE001` reading and failure-safe writing |
| Command-line workflow | File-to-file Haar denoising along the spectral axis |

## Wavelet coefficients

The coefficient catalogue contains filters for:

- Haar
- Beylkin
- Coiflet
- Daubechies
- Symmlet
- Vaidyanathan
- Battle-Lemarié

The normalized Haar coefficients are used by the validated one-level and multilevel Haar contracts. The other coefficient sets have not yet been integrated into a forward and inverse transform pipeline.

## Numerical contracts

Public dimensions use `size_t`.

Functions that allocate output return caller-owned memory. The caller is responsible for releasing returned allocations using the corresponding documented cleanup operation.

Invalid arguments and unsupported operations fail deterministically:

- Pointer-returning functions return `NULL`.
- Status-returning functions return `-1`.
- `errno` identifies errors such as `EINVAL`, `EOVERFLOW`, or `ENOTSUP`.

Operators do not modify caller-owned input matrices or filter coefficients unless explicitly documented.

For odd-length downsampling, samples at indices `0, 2, 4, ...` are retained. The output length is:

```text
ceil(input_length / 2)
```

Periodic convolution uses direct modular indexing and supports filters that are shorter than, equal to, or longer than the input signal.

## Repository structure

```text
.
├── Circular-Left-Shift/
├── Circular-Right-Shift/
├── Command-Line/
├── Discrete-Wavelet-Transform/
├── Haar-Denoising/
├── High-Pass-Downsampling-Operator/
├── High-Pass-Upsampling-Operator/
├── Low-Pass-Downsampling-Operator/
├── Low-pass-Upsampling-Operator/
├── Mirror-Filter/
├── Periodic-Convolution/
├── Rational-Transfer-Function/
├── Spectral-Cube/
├── Spectral-Cube-IO/
├── Time-Reversed-Periodic-Convolution/
├── Universal-Tools/
├── Upsampling-Operator/
├── Wavelet-Coefficients/
├── tests/
│   └── test_operators.c
├── .gitignore
├── main.c
├── makefile
└── README.md
```

`libspectral_operators.a` contains the numerical operators and spectral-cube file APIs. `tests/test_operators.c` verifies their public contracts, malformed-file handling, CLI validation, and deterministic failure behavior.

## Current limitations

The repository does not yet provide:

- Forward or inverse transforms for wavelets other than Haar
- A generic multidimensional transform contract
- Spatial-axis or full 3D cube transforms
- Automatic threshold selection
- Noise estimation
- Scientific formats other than `SCUBE001`, including ENVI, FITS, and HDF5
- Metadata handling
- Streaming and parallel processing
- End-to-end denoising quality validation

The CLI supports only Haar processing on the spectral axis with a user-provided threshold. The Haar paths have explicit pairing, normalization, level-limit, and per-level odd-length conventions. Phase, filter-ordering, boundary, axis, and subband-layout conventions remain undefined for generic and multidimensional transforms.

## Roadmap

1. Define and validate contracts for any additional wavelet families.
2. Add automatic threshold selection and noise estimation, then evaluate denoising quality.
3. Define spectral-cube metadata requirements before adding other scientific formats.
4. Stabilize the public API after the supported data and processing contracts are established.
5. Validate denoising behavior using documented datasets and quality metrics.

# License

Copyright (c) 2026 Jukka Joutsalainen

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE.
