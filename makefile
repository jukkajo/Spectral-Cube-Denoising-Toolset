CC = cc
AR = ar
CPPFLAGS = -I.
CFLAGS = -std=c11 -Wall -Wextra -Wpedantic -Wconversion -Wshadow
ARFLAGS = rcs
LDLIBS = -lm

TARGET = spectral-denoise
LIBRARY = libspectral_operators.a
TEST_TARGET = tests/test_operators
SANITIZER_TEST = tests/test_operators_sanitize
TEST_SOURCE = tests/test_operators.c
CLI_SOURCE = Command-Line/spectral_denoise_cli.c
CLI_OBJECT = $(CLI_SOURCE:.c=.o)

LIB_SOURCES = \
	Circular-Left-Shift/circular_left_shift.c \
	Circular-Right-Shift/circular_right_shift.c \
	Discrete-Wavelet-Transform/discrete_wavelet_transform.c \
	Haar-Denoising/haar_denoising.c \
	High-Pass-Downsampling-Operator/hi_pass_downsampling_operator.c \
	High-Pass-Upsampling-Operator/high_pass_upsampling_operator.c \
	Low-Pass-Downsampling-Operator/low_pass_downsampling_operator.c \
	Low-pass-Upsampling-Operator/low_pass_upsampling_operator.c \
	Mirror-Filter/mirror_filter.c \
	Periodic-Convolution/periodic_convolution.c \
	Rational-Transfer-Function/rational_transfer_function.c \
	Spectral-Cube/spectral_cube.c \
	Spectral-Cube-IO/spectral_cube_io.c \
	Time-Reversed-Periodic-Convolution/time_reversed_periodic_convolution.c \
	Universal-Tools/universal_tools.c \
	Upsampling-Operator/upsampling_operator.c \
	Wavelet-Coefficients/wavelet_coefficients.c

LIB_OBJECTS = $(LIB_SOURCES:.c=.o)
DEPENDENCIES = $(LIB_OBJECTS:.o=.d) $(CLI_OBJECT:.o=.d) main.d
SANITIZER_FLAGS = -fsanitize=address,undefined -fno-omit-frame-pointer
LEAK_CHECK ?= 0

.PHONY: all clean test sanitize

all: $(LIBRARY) $(TARGET)

$(TARGET): main.o $(CLI_OBJECT) $(LIBRARY)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ main.o $(CLI_OBJECT) $(LIBRARY) $(LDLIBS)

$(LIBRARY): $(LIB_OBJECTS)
	$(AR) $(ARFLAGS) $@ $^

$(TEST_TARGET): $(TEST_SOURCE) $(CLI_OBJECT) $(LIBRARY)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $(TEST_SOURCE) $(CLI_OBJECT) $(LIBRARY) $(LDLIBS)

test: $(TARGET) $(TEST_TARGET)
	./$(TEST_TARGET)

$(SANITIZER_TEST): $(TEST_SOURCE) $(CLI_SOURCE) $(LIB_SOURCES)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SANITIZER_FLAGS) $(LDFLAGS) -o $@ $(TEST_SOURCE) $(CLI_SOURCE) $(LIB_SOURCES) $(LDLIBS)

sanitize: $(TARGET) $(SANITIZER_TEST)
	# Use LEAK_CHECK=1 where LeakSanitizer is supported outside ptrace.
	ASAN_OPTIONS=detect_leaks=$(LEAK_CHECK):halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 ./$(SANITIZER_TEST)

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -MMD -MP -c -o $@ $<

clean:
	rm -f $(TARGET) $(LIBRARY) $(TEST_TARGET) $(SANITIZER_TEST) \
		program main.o $(CLI_OBJECT) $(LIB_OBJECTS) $(DEPENDENCIES)

-include $(DEPENDENCIES)
