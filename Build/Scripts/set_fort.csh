#!/bin/csh -f

# This script sets up the environment needed to compile Fortran routines.
# It expects either ia32 or intel64 as an input argument depengin on whether
# 32 or 64 bit routines are being compiled
#
# It assumes that an environment variable IFORT_COMPILER11 has been defined.
# This variable points to where the Intel Fortran compiler has been installed.
#
set platform=$1

if (! $?IFORT_COMPILER11) then
echo "***" Error: The environment variable IFORT_COMPILER11 is not defined.
echo Add an entry such as the following to one of your startup files
echo
echo setenv IFORT_COMPILER11 /opt/intel/Compiler/11.1/038
exit 1
endif
if (! -d $IFORT_COMPILER11) then
echo "*** Error: The Fortran compiler was not found.  It was expected at:"
echo "           $IFORT_COMPILER11 "
exit 1
endif

echo Using Fortran located at: $IFORT_COMPILER11
source $IFORT_COMPILER11/bin/ifortvars.csh $platform
set mpibin=
if ("$platform" == "ia32") then
set mpibin = /shared/openmpi_32/bin
endif
if ("$platform" == "intel64") then
if ($?INFINIBAND) then
set mpibin = /shared/openmpi_ib64/bin
else
set mpibin = /shared/openmpi_64/bin
endif
endif
echo "Using MPI distribution located at:$mpibin"
set path = ($mpibin $path)

exit 0
