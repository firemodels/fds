#!/bin/bash
IFORT_COMPILER /opt/intel/composerxe/

INTEL=$IFORT_COMPILER/bin/

source $INTEL/compilervars.sh intel64

make

