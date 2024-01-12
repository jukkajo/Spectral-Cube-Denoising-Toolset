CC = gcc
CFLAGS = -Wall -Wextra -std=c11
TARGET = program
SRC = main.c
LIB_DIR = ./Universal-Tools
LIB_SRC = $(LIB_DIR)/universal_tools.c

LIB_DIR2 = ./Discrete-Wavelet-Transform
LIB_SRC2 = $(LIB_DIR2)/discrete_wavelet_transform.c

INCLUDES = -I$(LIB_DIR) -I$(LIB_DIR2)

all: $(TARGET)

$(TARGET): $(SRC) $(LIB_SRC) $(LIB_SRC2)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^ -lm

clean:
	rm -f $(TARGET)

.PHONY: all clean

