#!/bin/bash
GEOMTEST=../intel_osx_64/geomtest
$GEOMTEST test.fds
$GEOMTEST test_azim.fds
$GEOMTEST test_elev.fds
$GEOMTEST test_scale.fds
