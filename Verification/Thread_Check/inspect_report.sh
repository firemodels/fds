#!/bin/bash

CURDIR=`pwd`
cd ../../..
GITROOT=`pwd`
cd $CURDIR

RESULT_DIR="inspect_results"
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

source /opt/intel19/inspector_2019/inspxe-vars.sh quiet



if [ "$showinput" == "1" ] ; then
  find . -name "$RESULT_DIR*"|while read fname; 
  do echo inspxe-cl -report $REPORT_TYPE -result-dir $fname
  done
  exit
fi

find . -name "$RESULT_DIR*"|while read fname;

do
echo $fname
 inspxe-cl -report $REPORT_TYPE -report-all -result-dir $fname;
sleep 3
done
