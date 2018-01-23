#!/bin/bash

if [ -e ~/FDS_SMV_ENVpc.sh ]; then
  source ~/FDS_SMV_ENVpc.sh
else
  source ~/FDS_SMV_ENV.sh
fi

export FDSEDITION=FDS6
export SMVEDITION=SMV6

export fds_smvroot=$linux_svn_root
export bundlebase=FDS_${fds_version}-SMV_${smv_version}_linux64
export runhost=$osx_hostname
export fdshost=$osx_hostname
export smvhost=$osx_hostname
export OSXBUNDLE=yes
export FDSVERSION=$fds_version
export SMVVERSION=$smv_version
export MPI_VERSION=$osx_mpi_version

export INSTALLDIR=FDS/$FDSEDITION

~/$fds_smvroot/fds/Build/Bundle/scripts/bundle_generic.sh
