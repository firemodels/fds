#!/bin/bash

source ../scripts/GET_ENV.sh

export fds_smvroot=$linux_svn_root
export bundlebase=${fds_version}-${smv_version}_osx64
export fdshost=$osx_hostname
export smvhost=$osx_hostname
export OSXBUNDLE=yes
export MPI_VERSION=$osx_mpi_version

export INSTALLDIR=FDS/$FDSEDITION

~/$fds_smvroot/fds/Build/Bundle/scripts/bundle_generic.sh
