#!/bin/bash

echo Building with Intel MPI
make -f ../makefile clean
make VPATH=".." -f ../makefile impi_intel_linux
