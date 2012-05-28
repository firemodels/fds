#/bin/bash -f

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

CURDIR=`pwd`
cd ..
export SVNROOT=`pwd`/..

# for OSX (without queuing)
#export BACKGROUND=$SVNROOT/Utilities/background/intel_osx_32/background
#export FDSEXE=$SVNROOT/FDS_Compilation/intel_osx_64/fds_intel_osx_64
#export FDS=$FDSEXE
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfds_noq.sh

# for Linux (with queing)
export FDSEXE=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDS=$FDSEXE

# blaze queue (default)
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds.sh

# fire60s queue
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfds6.sh
#export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds6.sh

# fire70s queue
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfds7.sh
#export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds7.sh

export BASEDIR=`pwd`
# uncomment following line to stop all cases
#export STOPFDS=1

scripts/SMV_Cases.sh

echo FDS cases submitted

