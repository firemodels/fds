#!/bin/bash
#
MAKEPO=../makepo/intel_linux_64/makepo_linux_64
SMVDIR=../source/smokeview
SHAREDDIR=../source/shared

echo updating smokeview_template.po
cat $SMVDIR/*.c $SMVDIR/*.cpp $SHAREDDIR/*.c | $MAKEPO | sort -u | $MAKEPO -c -a > smokeview_template.po

