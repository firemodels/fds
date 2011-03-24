#!/bin/bash -f
source $IFORT_COMIPLER/bin/iccvars.sh intel64
source $IFORT_COMIPLER/bin/ifortars.sh intel64
make -j4 -f ../Makefile intel_linux_64
