#!/bin/bash

# File name of the C source
SOURCE_FILE="Markov2D.c"
OUTPUT_FILE="markov_sim"

# Compilation command
echo "🔧 Compiling $SOURCE_FILE..."
gcc -O2 -Wall Markov2D.c -o markov_sim -lgsl -lgslcblas -lm
# Check if compilation succeeded
if [ $? -eq 0 ]; then
    echo "✅ Compilation successful! Run with:"
    echo "   ./markov_sim config.txt"
else
    echo "❌ Compilation failed."
fi
