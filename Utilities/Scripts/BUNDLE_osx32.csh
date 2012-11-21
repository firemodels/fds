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
setenv PLATFORM OSX32
setenv FDSEDITION $4
setenv FDSVERSION $5
setenv SMVVERSION $6
setenv MAJOR $7

setenv FDSOS _osx_32
setenv FDSOS2 _osx_32
setenv INSTALLDIR FDS/$FDSEDITION

$fds_smvroot/Utilities/Scripts/bundle_generic.csh
