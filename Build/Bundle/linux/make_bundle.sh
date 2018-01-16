#!/bin/bash
#
# this script is called from windows which passes in the directory 
# containing this script
#

export FDSEDITION=FDS6
export SMVEDITION=SMV6

export INTELLIBDIR=fire-notes/INSTALL/INTEL/INTEL_17u4/LIB
export INTELBINDIR=fire-notes/INSTALL/INTEL/INTEL_17u4/bin64
export OSLIBDIR=fire-notes/INSTALL/OSLIBS/LINUX

export fds_smvroot=$1
export bundlebase=$2
export runhost=$3
export fdshost=$3
export smvhost=$3
export FDSVERSION=$4
export SMVVERSION=$5
export MPI_VERSION=$6

export INSTALLDIR=FDS/$FDSEDITION
export MISCTO=LIB64
export COMPTO=INTEL

$fds_smvroot/fds/Build/Bundle/bundle_generic.sh
