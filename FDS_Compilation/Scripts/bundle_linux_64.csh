#!/bin/csh -f
set scriptdir=$1
set googledir=$scriptdir/../to_google
set makedir=$scriptdir/..

rm -f $googledir/fds5_intel_linux_64
rm -f $googledir/fds5_mpi_intel_linux_64
cp $makedir/intel_linux_64/fds5_intel_linux_64 $googledir/.
cp $makedir/mpi_intel_linux_64/fds5_mpi_intel_linux_64 $googledir/.
