#!/bin/bash
CURDIR=`pwd`
cd ..
export SVNROOT=`pwd`/..

export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export BASEDIR=`pwd`
# uncomment following line to stop all cases
#export STOPFDS=1


$RUNFDS Visualization2 bound_test
$RUNFDS Visualization2 cell_test
$RUNFDS Visualization2 diff2a
$RUNFDS Visualization2 diff2b
$RUNFDS Visualization2 fed_test2
$RUNFDS Visualization2 fed_test
$RUNFDS Visualization2 grid_test
$RUNFDS Visualization2 lava_lamp
$RUNFDS Visualization2 mplume5c2
$RUNFDS Visualization2 mplume5c2_large
$RUNFDS Visualization2 objects_dynamic2
$RUNFDS Visualization2 objects_dynamic3
$RUNFDS Visualization2 objects_dynamic4
$RUNFDS Visualization2 objects_include
$RUNFDS Visualization2 objects_movingbox
$RUNFDS Visualization2 objects_moving_slice
$RUNFDS Visualization2 objects_movingsphere2
$RUNFDS Visualization2 objects_rack
$RUNFDS Visualization2 objects_wui
$RUNFDS Visualization2 slice3d_test2
$RUNFDS Visualization2 slice3d_test
$RUNFDS Visualization2 smoke_test
$RUNFDS Visualization2 smv_reload
$RUNFDS Visualization2 vector_cell_test
