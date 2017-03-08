#!/bin/bash

./configure --prefix /shared/openmpi_64gnu CC=gcc CXX=g++ F77=gfortran FC=gfortran CFLAGS="-m64 -O2" CXXFLAGS="-m64 -O2" FFLAGS="-m64 -O2" FCFLAGS="-m64 -O2" LDFLAGS=-m64 --enable-static --disable-shared
