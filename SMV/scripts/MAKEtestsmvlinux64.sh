#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/intel_linux_64
./make_smv.sh -t
