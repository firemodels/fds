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
export OSXBUNDLE=yes
export PLATFORM=OSX64
export FDSEDITION=$4
export FDSVERSION=$5
export SMVVERSION=$6
export OPENMPI_VERSION=$7
export MAJOR=$8
export COMPLIBFROM=
export MISCLIBFROM=
export IB=

export FDSOS=_osx_64
export INSTALLDIR=FDS/$FDSEDITION

$fds_smvroot/fds/Utilities/Scripts/bundle_generic.sh
