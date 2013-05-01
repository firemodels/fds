#!/bin/bash

# inspect_report.sh
# Kristopher Overholt
# 5/1/2013

# Report results from thread checker

export SVNROOT=`pwd`/../..
source /opt/intel/inspector_xe/inspxe-vars.sh quiet

RESULT_DIR=$SVNROOT/Utilities/Scripts/inspect_openmp_ti3
REPORT_TYPE=problems

function usage {
echo "inspect_report.sh [ -d result-dir -r report-type ]"
echo "Report results from thread checker"
echo ""
echo "Options"
echo "-d result-dir - directory that contains thread checker results"
echo "-r report-type - type of report: problems [or] observations"
exit
}

while getopts 'd:hr:' OPTION
do
case $OPTION in
  d)
   cases="$OPTARG"
   ;;
  h)
   usage;
   ;;
  r)
   cases="$OPTARG"
   ;;
esac
done

inspxe-cl -report $REPORT_TYPE -result-dir $RESULT_DIR
