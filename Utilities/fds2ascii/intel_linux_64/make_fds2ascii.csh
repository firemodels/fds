#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh intel64

ifort -m64 -o fds2ascii_linux_64 ../../Data_Processing/fds2ascii.f90
