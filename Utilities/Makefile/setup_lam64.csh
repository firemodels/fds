#!/bin/csh -f
#
# script used to configure 64 bit version of LAM
#
# This script is run in the lam source direcory
#
setenv CFLAGS -m64
setenv CXXFLAGS -m64
setenv CC gcc
setenv CXX g++
setenv FC ifort
set DIR=/opt/intel/Compiler/11.0/074/bin/
#
# to configure a build for 32 bits change intel64 to ia32 in lnes below
#
source $DIR/iccvars.csh intel64
source $DIR/ifortvars.csh intel64
./configure --prefix=/usr/local/lam7i_64 -with-memory-manager=none --disable-shared
