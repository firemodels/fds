#!/bin/bash

# add -A to any case that you wish to be a part of the benchmark timing suite

$QFDS -d WUI BT10m_2x2km_LS.fds
#$QFDS -d WUI fire_line.fds
$QFDS -d WUI -A hill_structure.fds
$QFDS -d WUI levelset1.fds
#$QFDS -d WUI onetree_surf_1mesh.fds
$QFDS -d WUI pine_tree.fds
$QFDS -d WUI tree_test2.fds
$QFDS -d WUI wind_test1.fds
$QFDS -d WUI -A wind_test2.fds
