#!/bin/bash
source ../setopts.sh $*
make -f ../Makefile clean
make SMV_TESTFLAG="$SMV_TESTFLAG" SMV_TESTSTRING="$SMV_TESTSTRING" SMV_PROFILEFLAG="$SMV_PROFILEFLAG" SMV_PROFILESTRING="$SMV_PROFILESTRING" -f ../Makefile intel_osx_32
