#!/bin/bash

#  This script compiles a few utility programs.
 
source $IFORT_COMPILER/bin/compilervars.sh intel64

ifort -o fds2ascii fds2ascii.f90
ifort -o layer_height layer_height.f

