#!/bin/bash
#
# this script is called from windows which passes in the directory 
# containing this script
#

export fds_smvroot=$1
export bundlebase=$2
export runhost=$3
export fdshost=$3
export smvhost=$3
export PLATFORM=LINUX64
export FDSEDITION=$4
export FDSVERSION=$5
export SMVVERSION=$6
export MPI_VERSION=$7
export MAJOR=$8

export INTELLIBDIR=fire-notes/INSTALL/INTEL/INTEL_17u4/LIB
export INTELBINDIR=fire-notes/INSTALL/INTEL/INTEL_17u4/bin64
export OSLIBDIR=fire-notes/INSTALL/OSLIBS/LINUX

export FDSMODULE=$FDSEDITION

export FDSOS=_linux_64
export INSTALLDIR=FDS/$FDSEDITION
export MISCTO=LIB64
export COMPTO=INTEL

$fds_smvroot/fds/Utilities/Scripts/bundle_generic.sh
