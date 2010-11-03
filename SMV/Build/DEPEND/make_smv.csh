#!/bin/csh -f
make -f ../Makefile depend | grep -v gcc | grep -v g++ | grep -v gfortran > depend.d
