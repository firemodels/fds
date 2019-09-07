#!/bin/bash

echo Building intel_osx_64 version of fds2ascii
rm -f *.o
make -f ../Makefile intel_osx_64
