#!/bin/bash
CURDIR=`pwd`
cd ..
export SVNROOT=`pwd`/..

QFDS=/usr/local/bin/qfds.sh
# uncomment following line to stop all cases
#export STOPFDS=1

$QFDS -r -d Visualization2 bound_test.fds
$QFDS -r -d Visualization2 cell_test.fds
$QFDS -r -d Visualization2 diff2a.fds
$QFDS -r -d Visualization2 diff2b.fds
$QFDS -r -d Visualization2 fed_test2.fds
$QFDS -r -d Visualization2 fed_test.fds
$QFDS -r -d Visualization2 grid_test.fds
$QFDS -r -d Visualization2 lava_lamp.fds
$QFDS -r -d Visualization2 mplume5c2.fds
$QFDS -r -d Visualization2 mplume5c2_large.fds
$QFDS -r -d Visualization2 objects_dynamic2.fds
$QFDS -r -d Visualization2 objects_dynamic3.fds
$QFDS -r -d Visualization2 objects_dynamic4.fds
$QFDS -r -d Visualization2 objects_include.fds
$QFDS -r -d Visualization2 objects_movingbox.fds
$QFDS -r -d Visualization2 objects_moving_slice.fds
$QFDS -r -d Visualization2 objects_movingsphere2.fds
$QFDS -r -d Visualization2 objects_rack.fds
$QFDS -r -d Visualization2 objects_wui.fds
$QFDS -r -d Visualization2 slice3d_test2.fds
$QFDS -r -d Visualization2 slice3d_test.fds
$QFDS -r -d Visualization2 smoke_test.fds
$QFDS -r -d Visualization2 smv_reload.fds
$QFDS -r -d Visualization2 vector_cell_test.fds
