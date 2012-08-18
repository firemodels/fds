#!/bin/bash
SMV_TESTFLAG=
SMV_TESTSTRING=
while getopts 't' OPTION
do
case $OPTION in
  t)
   SMV_TESTFLAG="-D pp_BETA"
   SMV_TESTSTRING="test_"
  ;;
esac
done
export SMV_TESTFLAG
export SMV_TESTSTRING
