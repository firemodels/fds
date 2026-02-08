#!/bin/bash

echo Building with OpenMPI
make -f ../makefile clean
make VPATH=".." -f ../makefile ompi_gnu_linux
