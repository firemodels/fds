#!/bin/csh -f
source ../../scripts/set_cfort.csh intel64
make -j4 -f ../Makefile intel_linux_64
