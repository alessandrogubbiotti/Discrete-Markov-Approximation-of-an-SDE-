# Compiler and flags
CC      = gcc
CFLAGS  = -O2 -Wall -Iinclude
LDFLAGS = -lgsl -lgslcblas -lm

# Project structure
SRC_DIR   = src
INC_DIR   = include
BUILD_DIR = build
TARGET    = markov_sim

# Source files and object files
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Default target
all: $(TARGET)

# Link step
$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

# Compile step
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean build files
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# Run after building
run: all
	./$(TARGET)

