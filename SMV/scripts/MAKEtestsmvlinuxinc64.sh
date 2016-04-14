#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/smokeview/intel_linux_64
./make_smv_inc.sh -t
