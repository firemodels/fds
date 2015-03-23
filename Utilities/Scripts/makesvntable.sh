#!/bin/bash
# run from svnroot/Validation
CURDIR=`pwd`
cd ..
SVNROOT=~/FDS-SMV
maketable=$SVNROOT/Validation/Process_All_Output.sh
cd $CURDIR
grep VDIR $maketable | awk 'BEGIN { FS = "/" } ; { system("./makesvnentry.sh  " $2) }'  | sort -t '&' -k 2
