#!/bin/bash

source ~/.bashrc_setfort intel64
rm *.o
make -f ../Makefile intel_osx_64
