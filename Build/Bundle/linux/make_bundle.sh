#!/bin/bash

source ../scripts/GET_ENV.sh

export INTELLIBDIR=/var/local/bundle/INTEL/INTEL_18u2/LIB
export INTELBINDIR=/var/local/bundle/INTEL/INTEL_18u2/bin64
export OSLIBDIR=/var/local/bundle/OSLIBS/LINUX

export fds_smvroot=$linux_svn_root
export bundlebase=${fds_version}-${smv_version}_linux64
export fdshost=$linux_hostname
export smvhost=$linux_hostname
export MPI_VERSION=$linux_mpi_version

export INSTALLDIR=FDS/$FDSEDITION

~/$fds_smvroot/fds/Build/Bundle/scripts/bundle_generic.sh
