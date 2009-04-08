#!/bin/csh -f
set localrootdir=$1
set remoterootdir=$2
set googledir=$localrootdir/Utilities/Scripts/to_google
set makedir=$remoterootdir/Utilities/Makefile

rm -f $googledir/fds5_intel_osx_64
rm -f $googledir/fds5_mpi_intel_osx_64
scp devi1.nist.gov\:$makedir/Intel_OSX_64/fds5_intel_osx_64 $googledir/.
# scp devi1.nist.gov\:$makedir/Mpi_Intel_OSX_32/fds5_mpi_intel_osx_32 $googledir/.
