#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv runhost $3
setenv fdshost bluesky.cfr.nist.gov
setenv smvhost bluesky.cfr.nist.gov
setenv OSXBUNDLE

setenv manifest manifest_osx_32.html

setenv smokezipdir INTEL_OSX_32
setenv smokediffdir INTEL_OSX_32
setenv fds5dir intel_osx_32
setenv fds5mpidir mpi_intel_osx_32

setenv smokeview smv5_osx_32
setenv smokezip smokezip_osx_32
setenv smokediff smokediff_osx_32

setenv fds5 fds5_intel_osx_32
setenv fds5mpi fds5_mpi_intel_osx_32

setenv fds5out fds5_osx_32
setenv fds5mpiout fds5_mpi_osx_32

setenv fds2ascii fds2ascii_osx_32
setenv fds2asciidir intel_osx_32

setenv PLATFORM OSX32

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


