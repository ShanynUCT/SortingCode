#!/bin/bash

# Step 1: Source the environment
echo "Step 1: Sourcing the 'my_root_env' environment..."
source /Users/shanyn/miniconda/bin/activate my_root_env

# Check if sourcing was successful
if [ $? -ne 0 ]; then
    echo "Error: Environment sourcing failed. Exiting."
    exit 1
fi

# Step 2: Compile the C++ program
echo "Step 2: Compiling the C++ program..."
clear && g++ -std=c++0x -O3 sort2labr.C -o exe2 `root-config --cflags --libs` -lSpectrum

# Check if compilation was successful
if [ $? -ne 0 ]; then
    echo "Error: Compilation failed. Exiting."
    exit 1
fi

# Step 3: Execute the compiled program with arguments
echo "Step 3: Executing the compiled program..."
./exe2 ../../../Exp/2022/220615/analysis/angle/water/run41/R77.root source

# Check if execution was successful
if [ $? -ne 0 ]; then
    echo "Error: Program execution failed. Exiting."
    exit 1
fi

echo "Script completed successfully."
