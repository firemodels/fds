#!/bin/bash -f

export SVNROOT=`pwd`/../..
export BASEDIR=`pwd`
export INDIR=Current_Results
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh


source ~/.bashrc_fds intel64

# to override FDSMPI defined in above script, remove comment
# from  line below and define your own FDSMPI location
#export FDSMPI=override FDSMPI defined in above script

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Beyler_Hood_propane-1353-19-0
$RUNFDS $INDIR Beyler_Hood_propane-1353-19--10
$RUNFDS $INDIR Beyler_Hood_propane-1353-19-5
$RUNFDS $INDIR Beyler_Hood_propane-1825-19-0
$RUNFDS $INDIR Beyler_Hood_propane-1825-19--10
$RUNFDS $INDIR Beyler_Hood_propane-1825-19-5
$RUNFDS $INDIR Beyler_Hood_propane-2430-19--10
$RUNFDS $INDIR Beyler_Hood_propane-2430-19-5
$RUNFDS $INDIR Beyler_Hood_propane-3152-19-0
$RUNFDS $INDIR Beyler_Hood_propane-792-19-5
$RUNFDS $INDIR Beyler_Hood_propane-821-19-0
$RUNFDS $INDIR Beyler_Hood_propane-821-19--10

