#!/bin/bash

CURDIR=`pwd`

GITROOT=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  GITROOT=$FDSSMV
fi

RESULT_DIR=$GITROOT/FDS/Utilities/Scripts/inspect_openmp_ti3
REPORT_TYPE=problems
showinput=

function usage {
echo "inspect_report.sh [-d result-dir -h -r repository root -R report-type -v output command]"
echo "Report results from thread checker"
echo ""
echo "Options"
echo "-d result-dir - directory containing thread checker results"
echo "   [default: $RESULT_DIR]"
echo "-h - display this message"
echo "-r repository root - FDS repository root directory"
echo "   [default: $GITROOT]"
echo "-R report-type - type of report: problems or observations"
echo "   [default: $REPORT_TYPE]"
echo "-v - list command used to report thread checking results"
exit
}

while getopts 'd:hr:R:v' OPTION
do
case $OPTION in
  d)
   RESULT_DIR="$OPTARG"
   ;;
  h)
   usage;
   ;;
  r)
   GITROOT="$OPTARG"
   RESULT_DIR=$GITROOT/FDS/Utilities/Scripts/inspect_openmp_ti3
   ;;
  R)
   REPORT_TYPE="$OPTARG"
   ;;
  v)
   showinput=1
   ;;
esac
done

# Report results from thread checker

source /opt/intel/inspector_xe/inspxe-vars.sh quiet

if [ "$showinput" == "1" ] ; then
  echo inspxe-cl -report $REPORT_TYPE -result-dir $RESULT_DIR
  exit
fi
inspxe-cl -report $REPORT_TYPE -result-dir $RESULT_DIR
