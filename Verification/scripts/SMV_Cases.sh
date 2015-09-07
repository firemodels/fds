#!/bin/bash

$QFDS -d Visualization cell_test.fds
$RUNCFAST -d Visualization cfast_test.in
$QFDS -d Visualization colorbar.fds
$QFDS -d Visualization colorconv.fds
$QFDS -d Visualization fed_test.fds
$QFDS -d Visualization mplume5c8.fds
$QFDS -d Visualization objects_dynamic.fds
$QFDS -d Visualization objects_elem.fds
$QFDS -d Visualization objects_static.fds
$QFDS -d Visualization plume5c.fds
#$QFDS -p 8 -d Visualization mplume5c8.fds
$QFDS -d Visualization plume5cdelta.fds
$QFDS -d Visualization plumeiso.fds
$QFDS -d Visualization plume5c_bounddef.fds
$QFDS -d Visualization script_test.fds
$QFDS -d Visualization script_slice_test.fds
$QFDS -d Visualization sillytexture.fds
$QFDS -d Visualization slicemask.fds
$QFDS -d Visualization smoke_sensor.fds
$QFDS -d Visualization smoke_test.fds
$QFDS -d Visualization smoke_test2.fds
$QFDS -d Visualization sprinkler_many.fds
$QFDS -d Visualization thouse5.fds
$QFDS -d Visualization thouse5delta.fds
$QFDS -d Visualization transparency.fds
$QFDS -d Visualization vcirctest.fds
$QFDS -d Visualization vcirctest2.fds
$RUNTFDS -d Visualization version.fds

$QFDS -d WUI BT10m_2x2km_LS.fds
#$QFDS -d WUI fire_line.fds
$QFDS -d WUI hill_structure.fds
$QFDS -d WUI levelset1.fds
#$QFDS -d WUI onetree_surf_1mesh.fds
$QFDS -d WUI pine_tree.fds
$QFDS -d WUI tree_test2.fds
$QFDS -d WUI wind_test1.fds
