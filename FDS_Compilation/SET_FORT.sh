#!/bin/bash

if [ x$IFORT_COMPILER == "x" ]
then
  echo The environment variable \$IFORT_COMPILER is not defined or
  echo is not pointing to a valid location. Define \$IFORT_COMPILER
  echo in your startup .bashrc file.
  echo Example: 
  echo export IFORT_COMPILER=/opt/intel/composerxe-2011.1.107
  echo Script aborted.
  exit 
fi
platform=$1

if [ "$platform" == "ib64" ]
then
source $IFORT_COMPILER/bin/ifortvars.sh intel64
fi

if [ "$platform" == "intel64" ]
then
source $IFORT_COMPILER/bin/ifortvars.sh intel64
fi

if [ "$platform" == "ia32" ]
then
source $IFORT_COMPILER/bin/ifortvars.sh ia32
fi
