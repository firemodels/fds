#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv runhost $3
setenv smvhost bluesky.cfr.nist.gov
setenv fdshost bluesky.cfr.nist.gov
setenv OSXBUNDLE

setenv manifest manifest_osx_64.html

setenv smokezipdir INTEL_OSX_64
setenv smokediffdir INTEL_OSX_64
setenv fds5dir intel_osx_64
setenv fds5mpidir mpi_intel_osx_64

setenv smokeview smv5_osx_64
setenv smokezip smokezip_osx_64
setenv smokediff smokediff_osx_64

setenv fds5 fds5_intel_osx_64
setenv fds5mpi fds5_mpi_intel_osx_64

setenv fds5out fds5_osx_64
setenv fds5mpiout fds5_mpi_osx_64

setenv fds2ascii fds2ascii_osx_64
setenv fds2asciidir intel_osx_64

setenv PLATFORM OSX64

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


