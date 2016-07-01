#!/bin/bash
CURDIR=`pwd`

GITROOT=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  GITROOT=$FDSSMV
fi
RESULT_DIR=$GITROOT/FDS/Utilities/Scripts/inspect_openmp_ti3

function usage {
  echo "Usage: inspect_openmp.sh [-r repository root] [-v] casename.dfs"
  echo ""
  echo " -d result-dir - directory containing thread checker results"
  echo "    [default: $RESULT_DIR"
  echo " -h display this message"
  echo " -r repository root - FDS repository root directory"
  echo "    [default: $GITROOT]"
  echo " -v   - list command that will be used to thread check"
  echo "input_file - input file"
  echo ""
  exit
}

if [ $# -lt 1 ]
then
  usage
fi

showinput=
while getopts 'd:hr:v' OPTION
do
case $OPTION  in
  d)
   RESULT_DIR="$OPTARG"
   ;;
  h)
   usage;
   ;;
  r)
   GITROOT="$OPTARG"
   ;;
  v)
   showinput=1
   ;;
esac
done
shift $(($OPTIND-1))
case=$1

# Perform OpenMP thread checking (locate deadlocks and data races)

source /opt/intel/inspector_xe/inspxe-vars.sh quiet

TARGET=$GITROOT/FDS/Build/intel_linux_64_inspect

if [ "$showinput" == "" ] ; then
  if [ -d $RESULT_DIR ] ; then
    rm -r $RESULT_DIR
  fi
  cd $TARGET
  make -f ../makefile clean
  ./make_fds.sh
fi

export OMP_NUM_THREADS=2

if [ "$showinput" == "1" ] ; then
  echo inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $RESULT_DIR \
          -search-dir src=$GITROOT/FDS/Source \
          -- $GITROOT/FDS/Build/intel_linux_64_inspect/fds_intel_linux_64_inspect $case
  exit
fi
cd $CURDIR
inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $RESULT_DIR \
          -search-dir src=$GITROOT/FDS/Source \
          -- $GITROOT/FDS/Build/intel_linux_64_inspect/fds_intel_linux_64_inspect $case
