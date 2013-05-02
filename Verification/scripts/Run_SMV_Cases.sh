#/bin/bash -f

# This script runs the Smokeview Verification Cases on a 
# Linux machine with a batch queuing system

queue=
size=64
DEBUG=
OPENMP=

function usage {
echo "Run_SMV_Cases.sh [-d -h -o -p -q queue_name -s ]"
echo "Runs Smokeview verification suite"
echo ""
echo "Options"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-o - run OpenMP version of FDS"
echo "-p size - platform size"
echo "     default: 64"
echo "     other options: 32"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-s - stop FDS runs"
exit
}

CURDIR=`pwd`
cd ..
export SVNROOT=`pwd`/..

while getopts 'dhop:q:s' OPTION
do
case $OPTION in
  d)
   DEBUG=_db
   ;;
  h)
  usage;
  ;;
  o)
  OPENMP=openmp_
  ;;
  p)
   size="$OPTARG"
   ;;
  q)
   queue="$OPTARG"
   ;;
  s)
   stop_cases=true
   export STOPFDS=1
   ;;
esac
#shift
done

if [ "$size" != "32" ]; then
  size=64
fi
size=_$size

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
  PLATFORM2=osx_32
else
  PLATFORM=linux$size
  PLATFORM2=linux_32
fi


export BACKGROUND=$SVNROOT/Utilities/background/intel_$PLATFORM2/background
export FDSEXE=$SVNROOT/FDS_Compilation/${OPENMP}intel_$PLATFORM$DEBUG/fds_${OPENMP}intel_$PLATFORM$DEBUG
export FDS=$FDSEXE
export CFAST=~/cfast/CFAST/intel_$PLATFORM/cfast6_$PLATFORM

# Set queue to submit cases to

if [ "$queue" != "" ]; then
   queue="-q $queue"
fi

export BASEDIR=`pwd`

# Remove output files (unless stop option is used)
if [[ ! stop_cases ]] ; then
  echo "Removing FDS/CFAST output files"
  export RUNCFAST="$SVNROOT/Verification/scripts/Remove_CFAST_Files.sh"
  export RUNFDS="$SVNROOT/Verification/scripts/Remove_FDS_Files.sh"
  scripts/SMV_Cases.sh
  echo "FDS/CFAST output files removed"
fi

# run cases    

export RUNCFAST="$SVNROOT/Utilities/Scripts/runcfast.sh $queue"
export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh $queue"

if [ $RUN_OPENMP ]; then
  export RUNFDS="$SVNROOT/Utilities/Scripts/runfdsopenmp.sh $queue" 
fi

echo "" | $FDSEXE 2> $SVNROOT/Manuals/SMV_User_Guide/SCRIPT_FIGURES/fds.version

scripts/SMV_Cases.sh

echo FDS cases submitted

