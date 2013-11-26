#!/bin/bash
GEOMTEST=../intel_linux_64/geomtest
$GEOMTEST test.fds
$GEOMTEST test_azim.fds
$GEOMTEST test_elev.fds
$GEOMTEST test_scale.fds
$GEOMTEST test_group.fds
