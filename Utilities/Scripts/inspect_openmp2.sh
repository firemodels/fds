#!/bin/bash
CURDIR=`pwd`

GITROOT=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  GITROOT=$FDSSMV
fi

function usage {
  echo "Usage: inspect_openmp.sh [-r repository root] [-v] casename.dfs"
  echo ""
  echo " -h display this message"
  echo " -r repository root - name and location of repository where FDS is located"
  echo "    [default: $GITROOT]"
  echo " -v   - list command used to thread check"
  echo "input_file - input file"
  echo ""
  exit
}

if [ $# -lt 1 ]
then
  usage
fi

showinput=
while getopts 'hr:v' OPTION
do
case $OPTION  in
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

INSPECT_DIR=$GITROOT/Utilities/Scripts/inspect_openmp_ti3

# Perform OpenMP thread checking (locate deadlocks and data races)

source /opt/intel/inspector_xe/inspxe-vars.sh quiet

TARGET=$GITROOT/FDS_Compilation/intel_linux_64_inspect

if [ "$showinput" == "" ] ; then
  if [ -d $INSPECT_DIR ] ; then
    rm -r $INSPECT_DIR
  fi
  cd $TARGET
  make -f ../makefile clean
  ./make_fds.sh
fi

export OMP_NUM_THREADS=2

if [ "$showinput" == "1" ] ; then
  echo inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $INSPECT_DIR \
          -search-dir src=$GITROOT/FDS_Source \
          -- $GITROOT/FDS_Compilation/intel_linux_64_inspect/fds_intel_linux_64_inspect $case
  exit
fi
cd $CURDIR
inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $INSPECT_DIR \
          -search-dir src=$GITROOT/FDS_Source \
          -- $GITROOT/FDS_Compilation/intel_linux_64_inspect/fds_intel_linux_64_inspect $case
