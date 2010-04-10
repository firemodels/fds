#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh intel64

ifort -o ../fds2ascii_osx_64 ../../Data_Processing/fds2ascii.f90
