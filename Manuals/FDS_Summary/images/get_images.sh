#!/bin/bash
ROOT=$1

if [ "$ROOT" == "" ]; then
  ROOT=../../../..
fi

cp $ROOT/fds/Manuals/FDS_User_Guide/SCRIPT_FIGURES/*.png         user/.
cp $ROOT/fds/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/*.png verification/.
