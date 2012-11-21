#!/bin/csh -f
#
set here=`pwd`
set MAKEPO=../makepo/intel_linux_64/makepo_linux_64
set SMVDIR=../source/smokeview
set SHAREDDIR=../source/shared

echo updating smokeview_template.po
cat $SMVDIR/*.c $SMVDIR/*.cpp $SHAREDDIR/*.c $SHAREDDIR/*.cpp | $MAKEPO | sort -u | $MAKEPO -c -a > smokeview_template.po

