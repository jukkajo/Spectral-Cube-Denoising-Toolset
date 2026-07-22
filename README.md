# Spectral Cube Denoising Toolset

A C11 numerical operator library intended as the foundation for wavelet-based spectral-cube denoising.

The current implementation provides tested signal-processing primitives, multidimensional array utilities, an immutable wavelet coefficient catalogue, and a reproducible build and test workflow. Forward and inverse discrete wavelet transforms, denoising, spectral-cube file handling, and a user-facing command-line interface are not yet implemented.

## Status

The numerical operator layer is implemented and verified.

- The project builds without compiler warnings.
- All operator implementations link exactly once into a static library.
- The test suite contains 326 deterministic checks.
- AddressSanitizer and UndefinedBehaviorSanitizer checks pass.
- The discrete wavelet transform API currently returns `ENOTSUP` until its mathematical conventions are defined and validated.

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
program
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
| Discrete wavelet transform | Public boundary exists but currently returns `ENOTSUP` |

## Wavelet coefficients

The coefficient catalogue contains filters for:

- Haar
- Beylkin
- Coiflet
- Daubechies
- Symmlet
- Vaidyanathan
- Battle-Lemarié

The catalogue is available to the operator library, but the coefficient sets have not yet been integrated into a validated forward and inverse transform pipeline.

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
├── Discrete-Wavelet-Transform/
├── High-Pass-Downsampling-Operator/
├── High-Pass-Upsampling-Operator/
├── Low-Pass-Downsampling-Operator/
├── Low-pass-Upsampling-Operator/
├── Mirror-Filter/
├── Periodic-Convolution/
├── Rational-Transfer-Function/
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

`libspectral_operators.a` contains the numerical operator implementations. `tests/test_operators.c` verifies the public contracts and deterministic failure behavior.

The current `program` executable is a minimal integration entry point. It is not yet a spectral-cube denoising CLI.

## Current limitations

The repository does not yet provide:

- A validated forward discrete wavelet transform
- An inverse discrete wavelet transform
- Multilevel wavelet decomposition
- Coefficient thresholding or shrinkage
- Noise estimation
- Spectral-cube data structures
- Scientific file input or output
- Metadata handling
- A production command-line interface
- End-to-end denoising quality validation

The phase, filter-ordering, boundary, axis, and subband-layout conventions for the transform must be defined before the DWT is implemented.

## Roadmap

1. Define and validate the forward DWT contract against known reference vectors.
2. Implement inverse transformation and verify numerical round trips.
3. Add multilevel decomposition and a defined coefficient representation.
4. Implement and evaluate coefficient-domain denoising.
5. Add spectral-cube representation, processing axes, and file input/output.
6. Add a stable public API and command-line workflow.
7. Validate the complete pipeline using documented datasets and quality metrics.

# License

Copyright (c) 2026 Jukka Joutsalainen

Permission to use, copy, modify, and/or distribute this software for any purpose with or without fee is hereby granted, provided that the above copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE.
