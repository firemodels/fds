#!/bin/bash

# add -A to any case that you wish to be a part of the benchmark timing suite

$QFDS -d Visualization -t -A mplume5c8_bench.fds
$QFDS -d Visualization -t -A plume5c_bench.fds
$QFDS -d Visualization -t -A thouse5_bench.fds

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
$QFDS -d Visualization version.fds
$QFDS -d Visualization version2.fds
