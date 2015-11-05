#!/bin/bash
ssh blaze \( cd /home4/gforney/FDS-SMV/Utilities/wind2fds/intel_linux_64 \; rm *.o \; source $IFORT_COMPILER/bin/compilervars.sh intel64 \; make -f ../Makefile intel_linux_64 \)
