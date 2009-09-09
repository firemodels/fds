#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2

setenv INTELLIB INTEL_LIB64
setenv smokeview smv5_linux_64
setenv smokezip smokezip_linux
setenv fds5 fds5_intel_linux_64
setenv fds5dir intel_linux_64
setenv fds5mpi fds5_mpi_intel_linux_64
setenv fds5mpidir mpi_intel_linux_64
setenv fds2ascii fds2ascii_linux

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


