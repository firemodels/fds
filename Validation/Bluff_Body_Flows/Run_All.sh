#/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

#$QFDS $DEBUG $QUEUE -p  68 -d $INDIR geom_sphere_Re50_40D.fds
#$QFDS $DEBUG $QUEUE -p  68 -d $INDIR geom_sphere_Re100_40D.fds
#$QFDS $DEBUG $QUEUE -p  68 -d $INDIR geom_sphere_Re150_40D.fds
#$QFDS $DEBUG $QUEUE -p  68 -d $INDIR geom_sphere_Re200_40D.fds
#$QFDS $DEBUG $QUEUE -p  82 -d $INDIR geom_sphere_Re300_60D.fds

echo FDS cases submitted
