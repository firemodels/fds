#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv fdshost acrux.cfr.nist.gov

setenv fdsdir intel_linux_64
setenv fdsmpidir mpi_intel_linux_64

setenv fds fds_intel_linux_64
setenv fdsmpi fds_mpi_intel_linux_64

setenv fdsout fds_linux_64
setenv fdsmpiout fds_mpi_linux_64

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


