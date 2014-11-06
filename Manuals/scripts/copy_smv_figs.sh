#!/bin/bash

CURDIR=`pwd`
HOST=blaze.el.nist.gov
FDSREPO=FDS-SMV
DIRS="SMV_User_Guide SMV_Verification_Guide"

for dir in $DIRS
do
  scp $HOST\:~firebot/$FDSREPO/Manuals/$dir/SCRIPT_FIGURES/* ~/$FDSREPO/Manuals/$dir/SCRIPT_FIGURES/.
done

cd $CURDIR
