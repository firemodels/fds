#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh ia32

ifort -o ../fds2ascii_osx_32 ../../Data_Processing/fds2ascii.f90
