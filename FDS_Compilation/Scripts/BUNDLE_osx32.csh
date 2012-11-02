#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv runhost $3
setenv fdshost bluesky
setenv smvhost bluesky
setenv OSXBUNDLE

setenv manifest manifest_osx_32.html

setenv smokezipdir INTEL_OSX_32
setenv smokediffdir INTEL_OSX_32
setenv fdsdir intel_osx_32
setenv fdsmpidir mpi_intel_osx_32

setenv smokeview smv5_osx_32
setenv smokezip smokezip_osx_32
setenv smokediff smokediff_osx_32

setenv fds fds_intel_osx_32
setenv fdsmpi fds_mpi_intel_osx_32

setenv fdsout fds_osx_32
setenv fdsmpiout fds_mpi_osx_32

setenv fds2ascii fds2ascii_osx_32
setenv fds2asciidir intel_osx_32

setenv PLATFORM OSX32

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


