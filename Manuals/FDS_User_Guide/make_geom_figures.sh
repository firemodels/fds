#!/bin/bash
CURDIR=`pwd`

cd ../../Verification/scripts
./Run_SMV_Cases.sh -g

cd $CURDIR
cd ../../Verification/scripts
#call Make_SMV_Pictures.sh -g
cd $CURDIR

