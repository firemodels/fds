#!/bin/bash
../clean.sh $*
make -f ../Makefile intel_osx_test_64_dbg
