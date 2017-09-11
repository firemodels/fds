#!/bin/bash

CURDIR=`pwd`
cd ../../..
GITROOT=`pwd`
cd $CURDIR

RESULT_DIR=$CURDIR/inspect_results
REPORT_TYPE=problems
showinput=

function usage {
echo "inspect_report.sh -h -R report-type -v output command]"
echo "Report results from thread checker"
echo ""
echo "Options"
echo "-h - display this message"
echo "-R report-type - type of report: problems or observations"
echo "   [default: $REPORT_TYPE]"
echo "-v - list command used to report thread checking results"
exit
}

while getopts 'd:hr:R:v' OPTION
do
case $OPTION in
  h)
   usage;
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

source /opt/intel/inspector/inspxe-vars.sh quiet

if [ "$showinput" == "1" ] ; then
  echo inspxe-cl -report $REPORT_TYPE -result-dir $RESULT_DIR
  exit
fi
inspxe-cl -report $REPORT_TYPE -result-dir $RESULT_DIR
