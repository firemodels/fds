#!/bin/bash
ARG=$1

if [ "$ARG" == "" ]; then
# fortran type not specified on command line
# use INTEL_IFORT defined in .bashrc if it exists or
# ifort if it does not exist
   if [ "$INTEL_IFORT" == "" ]; then
      export INTEL_IFORT=ifort
   fi
else
# fortran type was specified on command line
# use it even if INTEL_IFORT was defined in .bashrc
# ifort if it does not
   export INTEL_IFORT=$ARG
fi
