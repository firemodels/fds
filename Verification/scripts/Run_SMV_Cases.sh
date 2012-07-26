#/bin/bash -f

# This script runs the Smokeview Verification Cases on a 
# Linux machine with a batch queuing system

queue=

function usage {
echo "Run_SMV_Cases.sh [-d -h -q queue_name -s ]"
echo "Runs Smokeview verification suite"
echo ""
echo "Options"
echo "-d - use debug version of FDS"
echo "-h - display this message"
echo "-q queue_name - run cases using the queue queue_name"
echo "     default: batch"
echo "     other options: fire60s, fire70s, vis"
echo "-s - stop FDS runs"
exit
}

CURDIR=`pwd`
cd ..
export SVNROOT=`pwd`/..

# for Linux (with queing)
# Set paths to FDS executable
# If no argument is specfied, then run FDS release version.

export FDSEXE=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDS=$FDSEXE
export CFAST=~/cfast/CFAST/intel_linux_64/cfast6_linux_64

# Otherwise, if -d (debug) option is specified, then run FDS DB version.
while getopts 'dhq:s' OPTION
do
case $OPTION in
  d)
   export FDSEXE=$SVNROOT/FDS_Compilation/intel_linux_64_db/fds_intel_linux_64_db
   export FDS=$FDSEXE
   ;;
  h)
  usage;
  ;;
  q)
   queue="$OPTARG"
   ;;
  s)
   export STOPFDS=1
   ;;
esac
shift
done

# Set queue to submit cases to

if [ "$queue" != "" ]; then
   queue="-q $queue"
fi
   
export RUNCFAST="$SVNROOT/Utilities/Scripts/runcfast.sh $queue"
export RUNFDS="$SVNROOT/Utilities/Scripts/runfds.sh $queue"
export RUNFDSFG="$SVNROOT/Utilities/Scripts/runfds.sh $queue"

export BASEDIR=`pwd`

scripts/SMV_Cases.sh

echo FDS cases submitted

