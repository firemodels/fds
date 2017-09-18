#!/bin/bash
#
# this script is called from windows which passes in the directory 
# containing this script
#
<<<<<<< HEAD
<<<<<<< HEAD
=======
INTEL_VERSION=17
=======
>>>>>>> cbc84ab0150a844a157a99ea0328c28921690a8d

>>>>>>> 900e38896dccc585db7f3099c98411e6fafb3e85
export fds_smvroot=$1
export bundlebase=$2
export runhost=$3
export fdshost=$3
export smvhost=$3
export PLATFORM=LINUX64
export FDSEDITION=$4
export FDSVERSION=$5
export SMVVERSION=$6
export OPENMPI_VERSION=$7
export MAJOR=$8
export COMPFROM=$9

export MISCFROM=fire-notes/INSTALL/LIBS/LINUX/LIB64
export IB=ib
export FDSMODULE=$FDSEDITION

export FDSOS=_linux_64
export INSTALLDIR=FDS/$FDSEDITION
export MISCTO=LIB64
<<<<<<< HEAD
<<<<<<< HEAD
export COMPTO=INTELLIBS16
=======
export COMPTO=INTELLIBS$INTEL_VERSION
>>>>>>> 900e38896dccc585db7f3099c98411e6fafb3e85
=======
export COMPTO=INTELLIBS
>>>>>>> cbc84ab0150a844a157a99ea0328c28921690a8d

$fds_smvroot/fds/Utilities/Scripts/bundle_generic.sh
