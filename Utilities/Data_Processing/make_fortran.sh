#!/bin/bash

#  This script compiles a few utility programs.
 
source $IFORT_COMPILER/bin/compilervars.sh intel64

echo use script in Utilities/fds2ascii/intel_linux_64 (or intel_osx_64 ) to build fds2ascii
ifort -o layer_height layer_height.f
ifort -o coupling_emulator coupling_emulator.f90

