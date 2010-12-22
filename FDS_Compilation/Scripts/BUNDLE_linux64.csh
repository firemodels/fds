#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv fdshost acrux.cfr.nist.gov

setenv fds5dir intel_linux_64
setenv fds5mpidir mpi_intel_linux_64

setenv fds5 fds5_intel_linux_64
setenv fds5mpi fds5_mpi_intel_linux_64

setenv fds5out fds5_linux_64
setenv fds5mpiout fds5_mpi_linux_64

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


