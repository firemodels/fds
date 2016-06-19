#!/bin/bash
#
# this script is called from windows which passes in the directory containing this script
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
export MAJOR=$7

export FDSOS=_linux_64
export INSTALLDIR=FDS/$FDSEDITION
export INTELLIB=~/FIRE-LOCAL/LIBS/LINUX/LIB64
export DESTLIB=LIB64

$fds_smvroot/Utilities/Scripts/bundle_generic.sh
