#!/bin/bash
# Run from base dir

export LD_LIBRARY_PATH=./include
make ./include/test_asymptotics.so && valgrind julia -i ./julia/interactive_test.jl

