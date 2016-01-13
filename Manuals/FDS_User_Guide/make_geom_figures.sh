#!/bin/bash
CURDIR=`pwd`

cd ../../Verification/scripts
./Run_SMV_Cases.sh -g  -w -j GM_

cd $CURDIR
cd ../../Verification/scripts
./Make_SMV_Pictures.sh -g
cd $CURDIR

