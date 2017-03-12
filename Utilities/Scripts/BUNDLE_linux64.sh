#!/bin/bash
#
# this script is called from windows which passes in the directory 
# containing this script
#
INTEL_VERSION=17

export fds_smvroot=$1
export bundlebase=$2
export runhost=$3
export fdshost=$3
export smvhost=$3
export PLATFORM=LINUX64
export FDSEDITION=$4
export FDSVERSION=$5
export SMVVERSION=$6
export MAJOR=$7
export COMPFROM=$8
export MISCFROM=$9

export FDSOS=_linux_64
export INSTALLDIR=FDS/$FDSEDITION
export MISCTO=LIB64
export COMPTO=INTELLIBS$INTEL_VERSION

$fds_smvroot/fds/Utilities/Scripts/bundle_generic.sh
