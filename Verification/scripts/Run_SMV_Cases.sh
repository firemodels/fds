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
# Set paths to FDS executable
# If no argument is specfied, then run FDS release version.
export FDSEXE=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export FDS=$FDSEXE

# Otherwise, if -d (debug) option is specified, then run FDS DB version.
while getopts 'd' OPTION
do
case $OPTION in
  d)
   export FDSEXE=$SVNROOT/FDS_Compilation/intel_linux_64_db/fds_intel_linux_64_db
   export FDS=$FDSEXE
   ;;
esac
done

# blaze queue (default)
export RUNCFAST=$SVNROOT/Utilities/Scripts/runcfast.sh
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds.sh

# fire60s queue
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfds6.sh
#export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds6.sh

# fire70s queue
#export RUNCFAST=$SVNROOT/Utilities/Scripts/runcfast7.sh
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfds7.sh
#export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds7.sh

# vis queue
#export RUNCFAST=$SVNROOT/Utilities/Scripts/runcfastv.sh
#export RUNFDS=$SVNROOT/Utilities/Scripts/runfdsv.sh
#export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfdsv.sh

export BASEDIR=`pwd`
# uncomment following line to stop all cases
#export STOPFDS=1

scripts/SMV_Cases.sh

echo FDS cases submitted

