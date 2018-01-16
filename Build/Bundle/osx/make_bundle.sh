#!/bin/bash
#
# this script is called from windows which passes in the directory 
# containing this script
#

export FDSEDITION=FDS6
export SMVEDITION=SMV6

export fds_smvroot=$1
export bundlebase=$2
export runhost=$3
export fdshost=$3
export smvhost=$3
export OSXBUNDLE=yes
export FDSVERSION=$4
export SMVVERSION=$5
export MPI_VERSION=$6
export COMPLIBFROM=
export MISCLIBFROM=

export INSTALLDIR=FDS/$FDSEDITION

$fds_smvroot/fds/Build/Bundle/scripts/bundle_generic.sh
