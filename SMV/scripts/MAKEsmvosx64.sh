#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/smokeview/intel_osx_64
./make_smv.sh
