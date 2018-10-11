#!/bin/bash

echo Building with Intel MPI
make -j4 VPATH=".." -f ../makefile impi_intel_linux
