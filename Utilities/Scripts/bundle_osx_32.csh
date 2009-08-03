#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot ~/$1
setenv scp_fds_smvroot $1
setenv makedir $scp_fds_smvroot/FDS_Compilation
setenv googledir $fds_smvroot/Utilities/to_google
setenv bundlebase $2
setenv bundledir $bundlebase/FDS/FDS5
setenv bundle_setup $fds_smvroot/Utilities/Scripts/bundle_setup
setenv mandir $fds_smvroot/Manuals/All_PDF_Files
setenv smvbindir $scp_fds_smvroot/SMV_5/bin
setenv fdshost tiger.cfr.nist.gov

setenv smokeview smv5_osx_32
setenv smokezip smokezip_osx
setenv fds5 fds5_intel_osx_32
setenv fds5dir intel_osx_32
setenv fds5mpi fds5_mpi_intel_osx_32
setenv fds5mpidir mpi_intel_osx_32
setenv fds2ascii fds2ascii_osx
setenv fds2asciidir $fds_smvroot/Utilities/fds2ascii

$fds_smvroot/Utilities/Scripts/bundle_generic.csh


