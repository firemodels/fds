#/bin/bash -f

# This script runs the Smokeview Verification Cases on a 
# Linux machine with a batch queuing system

queue=
size=64
DEBUG=

function usage {
echo "Run_SMV_Cases.sh [-d -h -p -q queue_name -s ]"
echo "Runs Smokeview verification suite"
echo ""
echo "Options"
echo "-d - use debug version of FDS"
echo "-h - display this message"
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

while getopts 'dhp:q:s' OPTION
do
case $OPTION in
  d)
   DEBUG=_db
   ;;
  h)
  usage;
  ;;
  p)
   size="$OPTARG"
   ;;
  q)
   queue="$OPTARG"
   ;;
  s)
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
export FDSEXE=$SVNROOT/FDS_Compilation/intel_$PLATFORM$DEBUG/fds_intel_$PLATFORM$DEBUG
export FDS=$FDSEXE
export CFAST=~/cfast/CFAST/intel_$PLATFORM/cfast6_$PLATFORM

# Set queue to submit cases to

if [ "$queue" != "" ]; then
   queue="-q $queue"
fi
   
export RUNCFAST="$SVNROOT/Utilities/Scripts/runcfast.sh $queue"
export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh $queue"

export BASEDIR=`pwd`

scripts/SMV_Cases.sh

echo FDS cases submitted

