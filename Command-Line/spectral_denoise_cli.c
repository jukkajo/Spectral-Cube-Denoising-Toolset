#include "spectral_denoise_cli.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../Spectral-Cube-IO/spectral_cube_io.h"

struct CliOptions {
    const char *input_path;
    const char *output_path;
    size_t levels;
    double threshold;
    enum HaarThresholdMode mode;
    int has_input;
    int has_output;
    int has_levels;
    int has_threshold;
    int has_mode;
};

static void print_usage(FILE *stream) {
    fprintf(
        stream,
        "Usage: spectral-denoise --input PATH --output PATH --levels N "
        "--threshold VALUE --mode hard|soft\n"
        "       spectral-denoise --help\n"
        "       spectral-denoise --version\n"
    );
}

static int argument_error(const char *message) {
    fprintf(stderr, "spectral-denoise: %s\n", message);
    return EXIT_FAILURE;
}

static int operation_error(const char *operation, const char *path) {
    const int saved_errno = errno;
    fprintf(
        stderr,
        "spectral-denoise: %s '%s': %s\n",
        operation,
        path,
        strerror(saved_errno)
    );
    errno = saved_errno;
    return EXIT_FAILURE;
}

static int parse_levels(const char *text, size_t *levels) {
    if (text == NULL || text[0] == '\0') {
        return -1;
    }

    size_t value = 0U;
    for (size_t index = 0U; text[index] != '\0'; ++index) {
        if (text[index] < '0' || text[index] > '9') {
            return -1;
        }
        const size_t digit = (size_t)(text[index] - '0');
        if (value > (SIZE_MAX - digit) / 10U) {
            return -1;
        }
        value = value * 10U + digit;
    }
    if (value == 0U) {
        return -1;
    }

    *levels = value;
    return 0;
}

static int parse_threshold(const char *text, double *threshold) {
    if (text == NULL || text[0] == '\0') {
        return -1;
    }

    char *end = NULL;
    errno = 0;
    const double value = strtod(text, &end);
    if (errno == ERANGE || end == text || *end != '\0' ||
        !isfinite(value) || value < 0.0) {
        return -1;
    }

    *threshold = value;
    return 0;
}

static int require_value(
    int argc,
    char *const argv[],
    int *index,
    const char **value
) {
    if (*index + 1 >= argc) {
        return -1;
    }
    ++(*index);
    *value = argv[*index];
    return 0;
}

static int parse_arguments(
    int argc,
    char *const argv[],
    struct CliOptions *options
) {
    for (int index = 1; index < argc; ++index) {
        const char *argument = argv[index];
        const char *value = NULL;

        if (strcmp(argument, "--input") == 0) {
            if (options->has_input != 0) {
                return argument_error("duplicate --input option");
            }
            if (require_value(argc, argv, &index, &value) != 0) {
                return argument_error("--input requires a path");
            }
            options->input_path = value;
            options->has_input = 1;
        } else if (strcmp(argument, "--output") == 0) {
            if (options->has_output != 0) {
                return argument_error("duplicate --output option");
            }
            if (require_value(argc, argv, &index, &value) != 0) {
                return argument_error("--output requires a path");
            }
            options->output_path = value;
            options->has_output = 1;
        } else if (strcmp(argument, "--levels") == 0) {
            if (options->has_levels != 0) {
                return argument_error("duplicate --levels option");
            }
            if (require_value(argc, argv, &index, &value) != 0 ||
                parse_levels(value, &options->levels) != 0) {
                return argument_error("--levels requires a positive integer");
            }
            options->has_levels = 1;
        } else if (strcmp(argument, "--threshold") == 0) {
            if (options->has_threshold != 0) {
                return argument_error("duplicate --threshold option");
            }
            if (require_value(argc, argv, &index, &value) != 0 ||
                parse_threshold(value, &options->threshold) != 0) {
                return argument_error(
                    "--threshold requires a finite nonnegative number"
                );
            }
            options->has_threshold = 1;
        } else if (strcmp(argument, "--mode") == 0) {
            if (options->has_mode != 0) {
                return argument_error("duplicate --mode option");
            }
            if (require_value(argc, argv, &index, &value) != 0) {
                return argument_error("--mode requires hard or soft");
            }
            if (strcmp(value, "hard") == 0) {
                options->mode = HAAR_THRESHOLD_HARD;
            } else if (strcmp(value, "soft") == 0) {
                options->mode = HAAR_THRESHOLD_SOFT;
            } else {
                return argument_error("--mode requires hard or soft");
            }
            options->has_mode = 1;
        } else if (strcmp(argument, "--help") == 0 ||
                   strcmp(argument, "--version") == 0) {
            return argument_error("--help and --version must be used alone");
        } else {
            return argument_error("unknown option");
        }
    }

    if (options->has_input == 0) {
        return argument_error("missing required --input option");
    }
    if (options->has_output == 0) {
        return argument_error("missing required --output option");
    }
    if (options->has_levels == 0) {
        return argument_error("missing required --levels option");
    }
    if (options->has_threshold == 0) {
        return argument_error("missing required --threshold option");
    }
    if (options->has_mode == 0) {
        return argument_error("missing required --mode option");
    }
    if (strcmp(options->input_path, options->output_path) == 0) {
        return argument_error("input and output paths must differ");
    }
    return EXIT_SUCCESS;
}

int spectral_denoise_cli_run(int argc, char *const argv[]) {
    if (argc == 2 && argv != NULL && strcmp(argv[1], "--help") == 0) {
        print_usage(stdout);
        return ferror(stdout) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (argc == 2 && argv != NULL && strcmp(argv[1], "--version") == 0) {
        puts(SPECTRAL_DENOISE_VERSION);
        return ferror(stdout) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    if (argc <= 0 || argv == NULL) {
        errno = EINVAL;
        return argument_error("invalid argument vector");
    }

    struct CliOptions options = {0};
    if (parse_arguments(argc, argv, &options) != EXIT_SUCCESS) {
        return EXIT_FAILURE;
    }

    SpectralCube *source = spectral_cube_read_file(options.input_path);
    if (source == NULL) {
        return operation_error("cannot read", options.input_path);
    }

    SpectralCube *destination = spectral_cube_denoise_spectral(
        source,
        options.levels,
        options.threshold,
        options.mode,
        NULL
    );
    if (destination == NULL) {
        const int saved_errno = errno;
        spectral_cube_free(source);
        errno = saved_errno;
        return operation_error("cannot denoise", options.input_path);
    }

    if (spectral_cube_write_file(options.output_path, destination) != 0) {
        const int saved_errno = errno;
        spectral_cube_free(destination);
        spectral_cube_free(source);
        errno = saved_errno;
        return operation_error("cannot write", options.output_path);
    }

    spectral_cube_free(destination);
    spectral_cube_free(source);
    return EXIT_SUCCESS;
}
