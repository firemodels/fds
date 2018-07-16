#!/bin/bash
CURDIR=`pwd`
cd ../../..
GITROOT=`pwd`
cd $CURDIR
RESULT_DIR=$CURDIR/inspect_results
OPENMPI=1

function usage {
  echo "Usage: inspect_openmp.sh [-v] casename.fds"
  echo ""
  echo " -I - use Intel MPI library"
  echo " -O - use Open MPI library [default]"
  echo " -h display this message"
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
while getopts 'hIOv' OPTION
do
case $OPTION  in
  h)
   usage;
   ;;
  I)
   OPENMPI=
   ;;
  O)
   OPENMPI=1
   ;;
  v)
   showinput=1
   ;;
esac
done
shift $(($OPTIND-1))
case=$1

# Perform OpenMP thread checking (locate deadlocks and data races)

source /opt/intel19/inspector_2019/inspxe-vars.sh quiet

if [ "$OPENMPI" == "1" ]; then
  INSPECTDIR=$GITROOT/fds/Build/mpi_intel_linux_64_inspect
  INSPECTAPP=$INSPECTDIR/fds_mpi_intel_linux_64_inspect
else
  INSPECTDIR=$GITROOT/fds/Build/impi_intel_linux_64_inspect
  INSPECTAPP=$INSPECTDIR/fds_impi_intel_linux_64_inspect
fi

if [ "$showinput" == "" ] ; then
  if [ -d $RESULT_DIR ] ; then
    rm -r $RESULT_DIR
  fi
  cd $INSPECTDIR
  make -f ../makefile clean
  ./make_fds.sh
fi

export OMP_NUM_THREADS=2

if [ "$showinput" == "1" ] ; then
  echo inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $RESULT_DIR \
          -search-dir src=$GITROOT/fds/Source \
          -- $INSPECTAPP $case
  exit
fi
cd $CURDIR
inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $RESULT_DIR \
          -search-dir src=$GITROOT/fds/Source \
          -- $INSPECTAPP $case
