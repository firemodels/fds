#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh ia32

ifort -m32 -o ../fds2ascii_linux_32 ../../Data_Processing/fds2ascii.f90
