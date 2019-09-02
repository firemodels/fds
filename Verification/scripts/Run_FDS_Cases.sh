#!/bin/bash

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

if [ ! -e .verification_script_dir ]; then
   echo "***error The Run_FDS_Cases.sh script must be run"
   echo "   in the directory Verification/script."
   echo "   script aborted."
   exit
fi

QUEUE=batch
DEBUG=
SINGLE=
nthreads=1
resource_manager=
walltime=
RUNOPTION=
CURDIR=`pwd`
QFDS_COUNT=/tmp/qfds_count_`whoami`
if [ "$BACKGROUND_PROG" == "" ]; then
  export BACKGROUND_PROG=background
fi
if [ "$BACKGROUND_DELAY" == "" ]; then
  BACKGROUND_DELAY=2
fi
if [ "$BACKGROUND_LOAD" == "" ]; then
  export BACKGROUND_LOAD=75
fi
REGULAR=1
BENCHMARK=1
OOPT=
POPT=
INTEL=i
INTEL2="-I"
GEOMCASES=1
INSPECTCASES=
WAIT=
EXE=
CHECKCASES=
RERUN=
DELAY=

function usage {
echo "Run_FDS_Cases.sh [ -d -h -m max_iterations -o nthreads -q queue_name "
echo "                   -s -r resource_manager ]"
echo "Runs FDS verification suite"
echo ""
echo "Options"
echo "-b - run only benchmark cases"
echo "-C - check that cases ran (used by firebot)"
echo "-d - use debug version of FDS"
echo "-D n - delay the submission of each case by n seconds"
echo "-e exe - run using exe"
echo "      Note: environment must be defined to use this executable"
echo "-F - rerun 'regular' cases that failed with 'BAD TERMINATION' errors"
echo "-g - run only geometry cases"
echo "-h - display this message"
echo "-j - job prefix"
echo "-J - use Intel MPI version of FDS"
echo "-m max_iterations - stop FDS runs after a specifed number of iterations (delayed stop)"
echo "     example: an option of 10 would cause FDS to stop after 10 iterations"
echo "-o nthreads - run FDS with a specified number of threads [default: $nthreads]"
echo "-O - use OpenMPI version of FDS"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire70s, vis"
echo "-r resource_manager - run cases using the resource manager"
echo "     default: PBS"
echo "     other options: SLURM"
echo "-R - run only regular (non-benchmark) cases"
echo "-s - stop FDS runs"
echo "-t - run only thread checking cases"
echo "-w time - walltime request for a batch job"
echo "     default: empty"
echo "     format for PBS: hh:mm:ss, format for SLURM: dd-hh:mm:ss"
echo "-W - wait for cases to complete before returning"
exit
}

function get_full_path {
  filepath=$1

  if [[ $filepath == /* ]]; then
    full_filepath=$filepath
  else
    dir_filepath=$(dirname  "${filepath}")
    filename_filepath=$(basename  "${filepath}")
    curdir=`pwd`
    cd $dir_filepath
    full_filepath=`pwd`/$filename_filepath
    cd $curdir
  fi
}

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
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} cases to complete."
        sleep 15
     done
   fi
}

cd ..
export SVNROOT=`pwd`/../..
cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR

while getopts 'bB:c:CdD:e:D:Fghj:JL:m:o:Oq:r:RsS:tw:W' OPTION
do
case $OPTION in
  b)
   BENCHMARK=1
   GEOMCASES=
   REGULAR=
   RERUN=
   ;;
  C)
   CHECKCASES="1"
   ;;
  d)
   DEBUG=_db
   SINGLE="1"
   ;;
  D)
   DELAY="-D $OPTARG"
   ;;
  e)
   EXE="$OPTARG"
   ;;
  F)
   BENCHMARK=
   GEOMCASES=
   REGULAR=
   RERUN=1
   ;;
  g)
   BENCHMARK=
   GEOMCASES=1
   REGULAR=
   RERUN=
   ;;
  h)
   usage;
   ;;
  j)
   JOBPREFIX="$OPTARG"
   ;;
  J)
   INTEL=i
   INTEL2="-I"
   ;;
  m)
   export STOPFDSMAXITER="$OPTARG"
   ;;
  o)
   nthreads="$OPTARG"
   ;;
  O)
   INTEL=
   INTEL2=
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   resource_manager="$OPTARG"
   ;;
  R)
   BENCHMARK=
   GEOMCASES=1
   REGULAR=1
   RERUN=
   ;;
  s)
   export STOPFDS=1
   ;;
  t)
   BENCHMARK=
   GEOMCASES=
   REGULAR=
   RERUN=
   INSPECTCASES=1
   DEBUG=_inspect
   ;;
  w)
   walltime="-w $OPTARG"
   ;;
  W)
   WAIT="1"
   ;;
esac
done

if [ "$JOBPREFIX" == "" ]; then
  JOBPREFIX=FB_
fi
export JOBPREFIX

size=_64

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
else
  PLATFORM=linux$size
fi
if [ "$OOPT" != "" ]; then
  OOPT="-O $OOPT"
fi
if [ "$POPT" != "" ]; then
  POPT="-O $POPT"
fi

if [ "$EXE" == "" ]; then
  export FDSMPI=$SVNROOT/fds/Build/${INTEL}mpi_intel_$PLATFORM$DEBUG/fds_${INTEL}mpi_intel_$PLATFORM$DEBUG
else
  get_full_path $EXE
  export FDSMPI=$full_filepath
fi

export QFDSSH="$SVNROOT/fds/Utilities/Scripts/qfds.sh $RUNOPTION $DELAY"

if [ "$resource_manager" == "SLURM" ]; then
   export RESOURCE_MANAGER="SLURM"
else
   export RESOURCE_MANAGER="PBS"
fi
if [ "$QUEUE" != "" ]; then
   if [ "$QUEUE" == "none" ]; then
      echo 0 > $QFDS_COUNT
   fi
   QUEUE="-q $QUEUE"
fi

export BASEDIR=`pwd`

export QFDS="$QFDSSH $walltime -n $nthreads $INTEL2 -e $FDSMPI $QUEUE $OOPT $POPT" 
if [ "$CHECKCASES" == "1" ]; then
  export QFDS="$SVNROOT/fds/Verification/scripts/Check_FDS_Cases.sh"
fi

cd ..
if [ "$BENCHMARK" == "1" ]; then
  if [ "$SINGLE" == "" ]; then
    ./FDS_Benchmark_Cases.sh
  else
    ./FDS_Benchmark_Cases_single.sh
  fi
  if [ "$CHECKCASES" == "" ]; then
    echo FDS benchmark cases submitted
  fi
fi

if [ "$CHECKCASES" == "1" ]; then
  export QFDS="$SVNROOT/fds/Verification/scripts/Check_FDS_Cases.sh"
else
  export QFDS="$QFDSSH $walltime -n $nthreads $INTEL2 -e $FDSMPI $QUEUE $OOPT $POPT" 
fi

cd $CURDIR
cd ..
if [ "$REGULAR" == "1" ]; then
  ./FDS_Cases.sh
  if [ "$CHECKCASES" == "" ]; then
    echo FDS non-benchmark cases submitted
  fi
fi
cd $CURDIR
cd ..
if [ "$GEOMCASES" == "1" ]; then
  ./GEOM_Cases.sh
  if [ "$CHECKCASES" == "" ]; then
    echo FDS geometry cases submitted
  fi
fi

cd $CURDIR
cd ..
if [ "$INSPECTCASES" == "1" ]; then
  ./INSPECT_Cases.sh
  if [ "$CHECKCASES" == "" ]; then
    echo FDS thread checking cases submitted
  fi
fi

cd $CURDIR
cd ..
if [ "$RERUN" == "1" ]; then
  grep 'BAD TERMINATION' */*.log | awk -F':' '{print($1)}' | sort -u | awk -F'/' '{print($2)}' | awk -F'.' '{print($1".fds")}' > badcaselist
  echo "#!/bin/bash" > RERUN_Cases.sh
  grep -f badcaselist FDS_Cases.sh >> RERUN_Cases.sh
  nlines=`cat RERUN_Cases.sh | wc -l`
  if [ $nlines -gt 1 ]; then
    echo warning the following cases failed with BAD TERMINATION errors. They were rerun
    grep 'BAD TERMINATION' -A 2 */*.log 
    chmod +x RERUN_Cases.sh
    ./RERUN_Cases.sh
    if [ "$CHECKCASES" == "" ]; then
      echo "FDS cases that failed with BAD TERMINATION errors re-submitted"
    fi
  fi
fi

if [ "$CHECKCASES" == "" ]; then
  if [ "$WAIT" == "1" ]; then
    wait_cases_end
  fi
fi
cd $CURDIR
