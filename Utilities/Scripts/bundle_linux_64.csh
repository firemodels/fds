#!/bin/csh -f
#
# this script is called from windows which passes in the directory containing this script
#
setenv fds_smvroot ~/$1
setenv makedir $fds_smvroot/FDS_Compilation
setenv googledir $fds_smvroot/Utilities/to_google
setenv bundlebase $2
setenv bundledir $bundlebase/FDS/FDS5
setenv bundle_setup $fds_smvroot/Utilities/Scripts/bundle_setup
setenv mandir $fds_smvroot/Manuals/All_PDF_Files
setenv smvbindir $fds_smvroot/SMV_5/bin
setenv smokeview smv5_linux_64
setenv smokezip smokezip_linux
setenv fds5 fds5_intel_linux_64
setenv fds5dir intel_linux_64
setenv fds5mpi fds5_mpi_intel_linux_64
setenv fds5mpidir mpi_intel_linux_64
setenv fds2ascii fds2ascii_linux
setenv fds2asciidir $fds_smvroot/Utilities/fds2ascii

$fds_smvroot/Utilities/Scripts/bundle_generic.csh

