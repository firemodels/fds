#!/bin/bash
./configure --prefix /shared/openmpi312_i19u4_64 \
  CC=icc CXX=icpc F77=ifort FC=ifort \
  CFLAGS="-m64 -O2" CXXFLAGS="-m64 -O2" \
  FFLAGS="-m64 -O2" FCFLAGS="-m64 -O2" \
  LDFLAGS=-m64 \
  --without-tm \
  --without-psm\
  --enable-mpirun-prefix-by-default \
  --without-verbs \
  --enable-static --disable-shared | tee CONFIGURE.out

make -j 6

