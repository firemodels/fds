#!/bin/bash
CURDIR=`pwd`
cd ../../..
GITROOT=`pwd`
cd $CURDIR
OPENMPI=

function usage {
  echo "Usage: build_inspect_openmp.sh [-v] casename.fds"
  echo ""
  echo " -I - use Intel MPI library [default]"
  echo " -O - use Open MPI library"
  echo " -h display this message"
  echo "input_file - input file"
  echo ""
  exit
}

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

cd $INSPECTDIR
make -f ../makefile clean
./make_fds.sh
