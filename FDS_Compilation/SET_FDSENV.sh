#!/bin/bash

# Copy SET_FDSENV.sh to SET_MYFDSENV.sh and edit the lines below beginning with FDS, FDSMPI
# FORTLIB and MPIDIST to match your distribution.

# This script is used to define the environment variables, FDS, FDSMPI, FORTLIB and MPIDIST,
# needed when running FDS. MPIDIST is also used when building a parallel FDS.  The script
# ifortvars sets up environment variables needed by the Fortran compiler.  
#
# This script is called by the various make_fds.sh scripts .  The script is called with one 
# argument, either ia32 (32 bit/ethernet), intel64 (64 bit/ethernet) or ib64 (64 bit/infiniband).
#
# This script is also called when running the VV scripts using the same input arguments.
#

platform=$1

if [ x$IFORT_COMPILER11 == "x" ]; then
echo "*** Error: The environment variable, IFORT_COMPILER11, is not defined.";
echo "    Define it in either ~/.cshrc or ~/.bashrc to point to the ";
echo "    Fortran compiler location.";
exit 1;
fi

if [ "$platform" == "ib64" ]
then
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64_ib/fds_intel_linux_64_ib
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_ib/fds_mpi_intel_linux_64_ib
export FORTLIB=/shared/LIB64
export MPIDIST=/shared/openmpi_ib64/
source $IFORT_COMPILER11/bin/ifortvars.sh intel64
fi

if [ "$platform" == "intel64" ]
then
export FDS=$SVNROOT/FDS_Compilation/mpi_intel_linux/fds_intel_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export FORTLIB=/shared/LIB64
export MPIDIST=/shared/openmpi_64/
source $IFORT_COMPILER11/bin/ifortvars.sh intel64
fi

if [ "$platform" == "ia32" ]
then
export FDS=$SVNROOT/FDS_Compilation/intel_linux_32/fds_intel_linux_32
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds_mpi_intel_linux_32
export FORTLIB=/shared/LIB32
export MPIDIST=/shared/openmpi_32/
source $IFORT_COMPILER11/bin/ifortvars.sh ia32
fi

export PATH=$MPIDIST/bin:$PATH
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh

echo ""
echo "*** FDS build/run environment ***"
echo "platform: $platform"
if [ x$SVNROOT != "x" ]; then
echo "serial FDS: $FDS"
echo "parallel FDS: $FDSMPI"
echo "runfdsmpi.sh: $RUNFDSMPI"
echo "Fortran libraries: $FORTLIB"
fi
echo "MPI distribution: $MPIDIST"
echo ""
