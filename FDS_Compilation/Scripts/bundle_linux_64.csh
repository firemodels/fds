#!/bin/csh -f
set scriptdir=$1
set googledir=$scriptdir/to_google
set makedir=$scriptdir/../Makefile

rm -f $googledir/fds5_intel_linux_64
rm -f $googledir/fds5_mpi_intel_linux_64
cp $makedir/Intel_Linux_64/fds5_intel_linux_64 $googledir/.
cp $makedir/Mpi_Intel_Linux_64/fds5_mpi_intel_linux_64 $googledir/.
