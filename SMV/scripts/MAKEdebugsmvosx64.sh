#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/intel_osx_64
./make_smv_db.sh
