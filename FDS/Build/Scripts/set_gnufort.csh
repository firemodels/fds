#!/bin/csh -f

# This script sets up the environment needed to compile Fortran routines.
# It expects either ia32 or intel64 as an input argument depengin on whether
# 32 or 64 bit routines are being compiled
#
# It assumes that an environment variable IFORT_COMPILER11 has been defined.
# This variable points to where the Intel Fortran compiler has been installed.
#
set platform=$1

echo Using GNU Fortran compilers
if ("$platform" == "ia32") then
set path = (/usr/local/lam7g/bin $path)
endif
if ("$platform" == "intel64") then
set path = (/usr/local/lam7g_64/bin $path)
endif

exit 0
