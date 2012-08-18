#!/bin/bash
SMV_TESTFLAG=
SMV_TESTSTRING=
SMV_PROFILEFLAG=
SMV_PROFILESTRING=
while getopts 'hpt' OPTION
do
case $OPTION in
  h)
  echo "options:"
  echo "-p - build a profiling version of smokeview"
  echo "-t - build a test version of smokeview"
  exit
  ;;
  p)
   SMV_PROFILEFLAG=" -p"
   SMV_PROFILESTRING="p"
  ;;
  t)
   SMV_TESTFLAG="-D pp_BETA"
   SMV_TESTSTRING="test_"
  ;;
esac
done
export SMV_TESTFLAG,SMV_TESTSTRING,SMV_PROFILEFLAG,SMV_PROFILESTRING
