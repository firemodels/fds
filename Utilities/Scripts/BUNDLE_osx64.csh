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
setenv PLATFORM OSX64
setenv FDSEDITION $4
setenv FDSVERSION $5
setenv SMVVERSION $6
setenv MAJOR $7

set OS=_osx_64

setenv manifest manifest$OS.html

setenv smokeviewdir intel$OS
setenv smokeview smokeview$OS
setenv smokeviewout smokeview$MAJOR$OS

setenv smokezipdir intel$OS
setenv smokezip smokezip$OS
setenv smokezipout smokezip$MAJOR$OS

setenv smokediffdir intel$OS
setenv smokediff smokediff$OS
setenv smokediffout smokediff$OS

setenv fdsdir intel$OS
setenv fds fds_intel$OS
setenv fdsout fds$MAJOR$OS

setenv fdsmpidir mpi_intel$OS
setenv fdsmpi fds_mpi_intel$OS
setenv fdsmpiout fds$MAJOR\_mpi$OS

setenv fds2asciidir intel$OS
setenv fds2ascii fds2ascii$OS
setenv fds2asciiout fds2ascii$MAJOR$OS

$fds_smvroot/Utilities/Scripts/bundle_generic.csh
