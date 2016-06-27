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
RUN_WUI=1
JOBPREFIX=
JOBPREF=
STOPFDS=
RUNOPTION=
CFASTREPO=~/cfastgitclean
COMPILER="intel"
WAIT=0

wait_cases_end()
{
   if [[ "$QUEUE" == "none" ]]
   then
     while [[ `ps -u $USER -f | fgrep .fds | grep -v grep` != '' ]]; do
        JOBS_REMAINING=`ps -u $USER -f | fgrep .fds | grep -v grep | wc -l`
        echo "Waiting for ${JOBS_REMAINING} cases to complete."
        sleep 15
     done
   else
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREF` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREF | wc -l`
        echo "Waiting for ${JOBS_REMAINING} cases to complete." 
        sleep 15
     done
   fi
}

function usage {
echo "Run_SMV_Cases.sh [-d -h -m max_iterations -o nthreads -p -q queue -s ]"
echo "Runs Smokeview verification suite"
echo ""
echo "Options"
echo "-c - cfast repo directory"
echo "-d - use debug version of FDS"
echo "-g - run only geometry cases"
echo "-h - display this message"
echo "-I - compiler (intel or gnu)"
echo "-j - job prefix"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-o nthreads - run OpenMP version of FDS with a specified number of threads [default: $nthreads]"
echo "-p size - platform size"
echo "     default: 64"
echo "     other options: 32"
echo "-q queue - run cases using the queue named queue"
echo "     default: batch"
echo "     other options: vis"
echo "-r - run only regular smokeview cases"
echo "-s - stop FDS runs"
echo "-u - use installed versions of utilities background and wind2fds"
echo "-w - wait for cases to complete before returning"
echo "-W - run only WUI cases"
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

export SVNROOT=`pwd`/../..
cd $SVNROOT
export SVNROOT=`pwd`

cd $CURDIR/..


use_installed="0"
while getopts 'c:dghI:j:m:o:p:q:rsuWw' OPTION
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
   RUN_GEOM=1
   RUN_WUI=0
   ;;
  h)
   usage;
   ;;
  I)
   COMPILER="$OPTARG"
   ;;
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  j)
   JOBPREFIX="-j $OPTARG"
   JOBPREF="$OPTARG"
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
   RUN_GEOM=0
   ;;
  s)
   stop_cases=true
   export STOPFDS=-s
   ;;
  u)
   use_installed="1"
   ;;
  w)
   WAIT="1"
   ;;
  W)
   RUN_SMV=0
   RUN_GEOM=0
   RUN_WUI=1
   ;;
esac
#shift
done

export FDS_DEBUG

size=_64

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
else
  PLATFORM=linux$size
fi
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

if [ "$use_installed" == "1" ] ; then
  export WIND2FDS=wind2fds
  export BACKGROUND=background
else
  export WIND2FDS=$SVNROOT/SMV/Build/wind2fds/${COMPILER}_$PLATFORM/wind2fds_$PLATFORM
  export BACKGROUND=$SVNROOT/SMV/Build/background/${COMPILER}_$PLATFORM/background
fi
export GEOM=$SVNROOT/SMV/source/geomtest/${COMPILER}_$PLATFORM/geomtest
export FDSEXE=$SVNROOT/FDS/Build/mpi_${COMPILER}_$PLATFORM$IB$DEBUG/fds_mpi_${COMPILER}_$PLATFORM$IB$DEBUG
export FDS=$FDSEXE
export FDSMPI=$SVNROOT/FDS/Build/mpi_${COMPILER}_$PLATFORM$IB$DEBUG/fds_mpi_${COMPILER}_$PLATFORM$IB$DEBUG
export CFAST=$CFASTREPO/Build/CFAST/${COMPILER}_$PLATFORM/cfast7_$PLATFORM
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
  scripts/GEOM_Cases.sh
  scripts/WUI_Cases.sh
  echo "FDS/CFAST output files removed"
fi

# run cases    

export  RUNCFAST="$QFDSSH -c -e $CFAST $QUEUE $STOPFDS $JOBPREFIX"
export      QFDS="$QFDSSH -e $FDSEXE $OPENMPOPTS $QUEUE $STOPFDS $JOBPREFIX"
export   RUNTFDS="$QFDSSH -e $FDSEXE $OPENMPOPTS $QUEUE $STOPFDS $JOBPREFIX"

echo "" | $FDSEXE 2> $SVNROOT/Manuals/SMV_User_Guide/SCRIPT_FIGURES/fds.version

if [[ ! $stop_cases ]] ; then
  if [ "$FDS_DEBUG" == "0" ] ; then
    if [ "$RUN_WUI" == "1" ] ; then
      is_file_installed $WIND2FDS
      cd $SVNROOT/Verification/WUI
      echo Converting wind data
      $WIND2FDS -prefix sd11 -offset " 100.0  100.0 0.0" wind_data1a.csv
    fi
  fi
fi

if [ "$RUN_SMV" == "1" ] ; then
  cd $SVNROOT/Verification
  scripts/SMV_Cases.sh
fi
if [ "$RUN_GEOM" == "1" ] ; then
  cd $SVNROOT/Verification
  scripts/GEOM_Cases.sh
fi
if [ "$RUN_WUI" == "1" ] ; then
  cd $SVNROOT/Verification
  scripts/WUI_Cases.sh
fi
if [ "$WAIT" == "1" ] ; then
  wait_cases_end
fi

echo FDS cases submitted

