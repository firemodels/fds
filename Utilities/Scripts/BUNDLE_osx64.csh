#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv runhost $3
setenv fdshost $3
setenv smvhost $3
setenv OSXBUNDLE
setenv PLATFORM OSX64
setenv FDSEDITION $4
setenv FDSVERSION $5
setenv SMVVERSION $6
setenv MAJOR $7

setenv FDSOS _osx_64
setenv FDSOS2 _osx_32
setenv INSTALLDIR FDS/$FDSEDITION

$fds_smvroot/Utilities/Scripts/bundle_generic.csh
