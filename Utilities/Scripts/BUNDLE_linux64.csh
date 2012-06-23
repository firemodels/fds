#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot $1
setenv bundlebase $2
setenv runhost $3
setenv fdshost $3
setenv smvhost $3
setenv FDSEDITION $4
setenv FDSVERSION $5
setenv FDSVERSION $6

setenv PLATFORM LINUX64

setenv manifest manifest_linux_64.html

setenv INTELLIB /shared/LIB64
setenv smokezipdir intel_linux_64
setenv smokediffdir intel_linux_64
setenv fdsdir intel_linux_64
setenv fdsmpidir mpi_intel_linux_64

setenv smokeviewdir intel_linux_64
setenv smokeview smokeview_linux_64
setenv smokezip smokezip_linux_64
setenv smokediff smokediff_linux_64

setenv fds fds_intel_linux_64
setenv fdsmpi fds_mpi_intel_linux_64

setenv fdsout fds_linux_64
setenv fdsmpiout fds_mpi_linux_64

setenv fds2ascii fds2ascii_linux_64
setenv fds2asciidir intel_linux_64

csh $fds_smvroot/Utilities/Scripts/bundle_generic.csh


