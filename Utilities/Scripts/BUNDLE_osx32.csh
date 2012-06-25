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

set OS=_osx_32

$fds_smvroot/Utilities/Scripts/bundle_generic.csh
