#!/bin/bash

source ../scripts/GET_ENV.sh

export FDSEDITION=FDS6
export SMVEDITION=SMV6

export INTELLIBDIR=/var/local/bundle/INTEL/INTEL_17u4/LIB
export INTELBINDIR=/var/local/bundle/INTEL/INTEL_17u4/bin64
export OSLIBDIR=/var/local/bundle/OSLIBS/LINUX

export fds_smvroot=$linux_svn_root
export bundlebase=FDS_${fds_version}-SMV_${smv_version}_linux64
export runhost=$linux_hostname
export fdshost=$linux_hostname
export smvhost=$linux_hostname
export FDSVERSION=$fds_version
export SMVVERSION=$smv_version
export MPI_VERSION=$linux_mpi_version

export INSTALLDIR=FDS/$FDSEDITION
export MISCTO=LIB64
export COMPTO=INTEL

~/$fds_smvroot/fds/Build/Bundle/scripts/bundle_generic.sh
