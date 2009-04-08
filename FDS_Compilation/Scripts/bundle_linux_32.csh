#!/bin/csh -f
set scriptdir=$1
set googledir=$scriptdir/to_google
set makedir=$scriptdir/../Makefile

rm -f $googledir/fds5_intel_linux_32
rm -f $googledir/fds5_mpi_intel_linux_32
cp -v $makedir/Intel_Linux_32/fds5_intel_linux_32 $googledir/.
cp -v $makedir/Mpi_Intel_Linux_32/fds5_mpi_intel_linux_32 $googledir/.
