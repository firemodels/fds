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
export OPENMPI_VERSION=$7
export MAJOR=$8
export COMPFROM=$9

#export MISCFROM=$10
export MISCFROM=fire-notes/INSTALL/LIBS/LINUX/LIB64
#export IB=$11 
export IB=ib

export FDSOS=_linux_64
export INSTALLDIR=FDS/$FDSEDITION
export MISCTO=LIB64
export COMPTO=INTELLIBS

$fds_smvroot/fds/Utilities/Scripts/bundle_generic.sh
