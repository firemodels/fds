#!/bin/bash
IFORT_COMPILER=/opt/intel/composerxe/

INTEL=$IFORT_COMPILER/bin/

source $INTEL/compilervars.sh intel64

./configure --prefix /shared/openmpi_64ib \
  CC=icc CXX=icpc F77=ifort FC=ifort CFLAGS="-m64 -O2" CXXFLAGS="-m64 -O2" FFLAGS="-m64 -O2" FCFLAGS="-m64 -O2" LDFLAGS=-m64 \
  --with-tm=/usr/local/torque \
  --enable-mpirun-prefix-by-default \
  --with-verbs=/usr \
  --with-verbs-libdir=/usr/lib64 --enable-static --disable-shared | tee CONFIGURE.out



