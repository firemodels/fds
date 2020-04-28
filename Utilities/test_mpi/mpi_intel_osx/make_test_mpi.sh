#!/bin/bash

echo Building with OpenMPI
make -f ../makefile clean
make VPATH=".." -f ../makefile mpi_intel_osx
