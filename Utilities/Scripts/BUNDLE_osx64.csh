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
setenv PLATFORM OSX64
setenv FDS_EDITION $4

setenv manifest manifest_osx_64.html

setenv smokezipdir intel_osx_64
setenv smokediffdir intel_osx_64
setenv fdsdir intel_osx_64
setenv fdsmpidir mpi_intel_osx_64

setenv smokeviewdir intel_osx_64
setenv smokeview smokeview_osx_64
setenv smokezip smokezip_osx_64
setenv smokediff smokediff_osx_64

setenv fds fds_intel_osx_64
setenv fdsmpi fds_mpi_intel_osx_64

setenv fdsout fds_osx_64
setenv fdsmpiout fds_mpi_osx_64

setenv fds2ascii fds2ascii_osx_64
setenv fds2asciidir intel_osx_64

$fds_smvroot/Utilities/Scripts/bundle_generic.csh
