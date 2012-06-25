#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv runhost $3
setenv fdshost blaze.nist.gov
setenv smvhost blaze.nist.gov
setenv PLATFORM LINUX64
setenv FDSEDITION $4
setenv FDSVERSION $5
setenv SMVVERSION $6
setenv MAJOR $7

setenv FDSOS _linux_64

$fds_smvroot/Utilities/Scripts/bundle_generic.csh
