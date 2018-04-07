#!/bin/bash
# Run from base dir

export LD_LIBRARY_PATH=../../include
make test_asymptotics.so && julia -i ./interactive_test.jl

