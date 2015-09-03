#!/bin/bash
# run from svnroot/Validation
CURDIR=`pwd`
cd ..
if [ "$FDSSMV" == "" ] ; then
  FDSSMV=~/FDS-SMVgitclean
fi
maketable=$FDSSMV/Validation/Process_All_Output.sh
cd $CURDIR
grep VDIR $maketable | awk 'BEGIN { FS = "/" } ; { system("./makegitentry.sh  " $2) }'  | sort -t '&' -k 2
