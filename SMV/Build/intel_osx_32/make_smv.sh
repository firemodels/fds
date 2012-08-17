#!/bin/bash
source ../clean.sh $*
make -f ../Makefile clean
make SMV_TESTFLAG="$SMV_TESTFLAG" SMV_TESTSTRING="$SMV_TESTSTRING" -f ../Makefile intel_osx_32
