#/bin/bash

# This script runs the Smokeview Verification Cases on a 
# Linux machine with a batch queuing system

QUEUE=batch
size=64
DEBUG=
OPENMP_OPTS=
FDS_DEBUG=0
nthreads=1
RUN_SMV=1
RUN_GEOM=1
JOBPREFIX=
# not running any mpi cases now
RUN_MPI=0
STOPFDS=
RUNOPTION=
CFASTREPO=cfastgitclean

function usage {
echo "Run_SMV_Cases.sh [-d -h -m max_iterations -o nthreads -p -q queue -s ]"
echo "Runs Smokeview verification suite"
echo ""
echo "Options"
echo "-c - cfast repo directory"
echo "-d - use debug version of FDS"
echo "-g - run only geometry cases"
echo "-h - display this message"
echo "-j - job prefix"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-M - run only cases using multiple processes"
echo "-o nthreads - run OpenMP version of FDS with a specified number of threads [default: $nthreads]"
echo "-p size - platform size"
echo "     default: 64"
echo "     other options: 32"
echo "-q queue - run cases using the queue named queue"
echo "     default: batch"
echo "     other options: vis"
echo "-r - run only regular smokeview cases"
echo "-s - stop FDS runs"
echo "-S - run only cases using a single process"
echo "-u - use installed versions of utilities background and wind2fds"
exit
}

is_file_installed()
{
  program=$1
  notfound=`$program -help 2>&1 | tail -1 | grep "not found" | wc -l`
  if [ "$notfound" == "1" ] ; then
    echo "***error: $program not available. Run aborted." 
    exit
  fi
}


CURDIR=`pwd`
cd ..

export SVNROOT=`pwd`/..
cd $SVNROOT
export SVNROOT=`pwd`

cd $CURDIR/..


use_installed="0"
while getopts 'c:dghj:Mm:o:p:q:rSsu' OPTION
do
case $OPTION in
  c)
   CFASTREPO="$OPTARG"
   ;;
  d)
   DEBUG=_db
   FDS_DEBUG=1
   ;;
  g)
   RUN_SMV=0
   RUN_MPI=0
   RUN_GEOM=1
   ;;
  h)
   usage;
   ;;
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  M)
   RUNOPTION="-M"
   ;;
  j)
   JOBPREFIX="-j $OPTARG"
   ;;
  o)
   nthreads="$OPTARG"
   OPENMP_OPTS="-n $nthreads"
   ;;
  p)
   size="$OPTARG"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   RUN_SMV=1
   RUN_MPI=0
   RUN_GEOM=0
   ;;
  s)
   stop_cases=true
   export STOPFDS=-s
   ;;
  S)
   RUNOPTION="-S"
   ;;
  u)
   use_installed="1"
   ;;
esac
#shift
done

export FDS_DEBUG

size=_64

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
  PLATFORM2=osx_32
  PLATFORM3=osx_64
else
  PLATFORM=linux$size
  PLATFORM2=linux_32
  PLATFORM3=linux_64
fi
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

if [ "$use_installed" == "1" ] ; then
  export WIND2FDS=wind2fds
  export BACKGROUND=background
else
  export WIND2FDS=$SVNROOT/Utilities/wind2fds/intel_$PLATFORM/wind2fds_$PLATFORM
  export BACKGROUND=$SVNROOT/Utilities/background/intel_$PLATFORM2/background
fi
export GEOM=$SVNROOT/SMV/source/geomtest/intel_$PLATFORM/geomtest
export FDSEXE=$SVNROOT/FDS_Compilation/mpi_intel_$PLATFORM$IB$DEBUG/fds_mpi_intel_$PLATFORM$IB$DEBUG
export FDS=$FDSEXE
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_$PLATFORM$IB$DEBUG/fds_mpi_intel_$PLATFORM$IB$DEBUG
export CFAST=~/$CFASTREPO/CFAST/intel_$PLATFORM/cfast7_$PLATFORM
QFDSSH="$SVNROOT/Utilities/Scripts/qfds.sh $RUNOPTION"

# Set queue to submit cases to

if [ "$QUEUE" != "" ]; then
   if [ "$QUEUE" == "none" ]; then
      is_file_installed $BACKGROUND
   fi
   QUEUE="-q $QUEUE"
fi

export BASEDIR=`pwd`

# Remove output files (unless stop option is used)
if [[ ! $stop_cases ]] ; then
  echo "Removing FDS/CFAST output files"
  export RUNCFAST="$SVNROOT/Verification/scripts/Remove_CFAST_Files.sh"
  export QFDS="$SVNROOT/Verification/scripts/Remove_FDS_Files.sh"
  export RUNTFDS="$SVNROOT/Verification/scripts/Remove_FDS_Files.sh"
  scripts/SMV_Cases.sh
  echo "FDS/CFAST output files removed"
fi

# run cases    

export  RUNCFAST="$QFDSSH -c -e $CFAST $QUEUE $STOPFDS $JOBPREFIX"
export      QFDS="$QFDSSH -e $FDSEXE $OPENMPOPTS $QUEUE $STOPFDS $JOBPREFIX"
export   RUNTFDS="$QFDSSH -e $FDSEXE $OPENMPOPTS $QUEUE $STOPFDS $JOBPREFIX"

echo "" | $FDSEXE 2> $SVNROOT/Manuals/SMV_User_Guide/SCRIPT_FIGURES/fds.version

if [[ ! $stop_cases ]] ; then
  if [ "$FDS_DEBUG" == "0" ] ; then

    is_file_installed $WIND2FDS
    cd $SVNROOT/Verification/WUI
    echo Converting wind data
    $WIND2FDS -prefix sd11 -offset " 100.0  100.0 0.0" wind_data1a.csv
  fi
fi

if [ "$RUN_SMV" == "1" ] ; then
  cd $SVNROOT/Verification
  scripts/SMV_Cases.sh
fi
if [ "$RUN_GEOM" == "1" ] ; then
  cd $SVNROOT/Verification
  scripts/SMV_geom_Cases.sh
fi
if [ "$RUN_MPI" == "1" ] ; then
  cd $SVNROOT/Verification
  export QFDS="$QFDSSH -e $FDSMPI $QUEUE $STOPFDS $JOBPREFIX"
  scripts/SMV_MPI_Cases.sh
fi

echo FDS cases submitted

