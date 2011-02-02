#!/bin/bash

platform=$1
option=$2

# To build FDS and/or run the various VV scripts do the following:
# 1.  Copy SET_FDSENV.sh to SET_MYFDSENV.sh
# 2.  In SET_MYFDSENV.sh edit lines below defining locations of
#     Fortran compiler and various MPI distributions.
# Note: you only need to define those MPI variables (if any) for
#       which you wish to build FDS versions for.

export IFORT_COMPILER11=/opt/intel/composerxe-2011.1.107/
export MPIDIST32=/shared/openmpi_32
export MPIDIST64=/shared/openmpi_64
export MPIDIST64IB=/shared/openmpi_64ib

#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# shouldn't need to edit any lines below

# --------- define environment for building FDS ------------------

if [ "$option" == "build" ]
then
if [ "$platform" == "ib64" ]
then
source $IFORT_COMPILER11/bin/ifortvars.sh intel64
export MPIDIST=$MPIDIST64IB
fi

if [ "$platform" == "intel64" ]
then
source $IFORT_COMPILER11/bin/ifortvars.sh intel64
export MPIDIST=$MPIDIST64
fi

if [ "$platform" == "ia32" ]
then
source $IFORT_COMPILER11/bin/ifortvars.sh ia32
export MPIDIST=$MPIDIST32
fi
fi

# ------------ define environment for running VV scripts FDS ---------

if [ "$option" == "run" ]
if [ "$platform" == "ib64" ]
then
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64_ib/fds_intel_linux_64_ib
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_ib/fds_mpi_intel_linux_64_ib
export MPIDIST=$MPIDIST64IB
fi

if [ "$platform" == "intel64" ]
then
export FDS=$SVNROOT/FDS_Compilation/mpi_intel_linux/fds_intel_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export MPIDIST=$MPIDIST64
fi

if [ "$platform" == "ia32" ]
then
export FDS=$SVNROOT/FDS_Compilation/intel_linux_32/fds_intel_linux_32
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds_mpi_intel_linux_32
export MPIDIST=$MPIDIST32
fi
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
fi

export PATH=$MPIDIST/bin:$PATH

if [ "$option" == "build" ]
echo ""
echo "*** FDS build environment ***"
echo "platform: $platform"
echo "MPI distribution: $MPIDIST"
echo ""
fi
if [ "$option" == "run" ]
echo ""
echo "*** FDS build environment ***"
echo "platform: $platform"
echo "serial FDS: $FDS"
echo "parallel FDS: $FDSMPI"
echo "MPI distribution: $MPIDIST"
echo ""
fi
