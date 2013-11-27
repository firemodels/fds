#!/bin/bash
GEOMTEST=../intel_osx_64/geomtest
$GEOMTEST test_azim.fds
$GEOMTEST test_component.fds
$GEOMTEST test_elev.fds
$GEOMTEST test_group.fds
$GEOMTEST test_group2.fds
$GEOMTEST test_group3.fds
$GEOMTEST test_scale.fds
