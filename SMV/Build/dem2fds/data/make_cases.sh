#!/bin/bash
option=$1 $2 
if [ "$option" == ""]; then
  option=-obst
fi

dem2fds=dem2fds
#dem2fds=../intel_linux_64/dem2fds_linux_64

terraindir=~/terrain

$dem2fds $option -dir $terraindir/blodget blodget.in 

$dem2fds $option -dir $terraindir/demtest demtest1.in 
$dem2fds $option -dir $terraindir/demtest demtest2.in 

$dem2fds $option -dir $terraindir/nist nist.in 

$dem2fds $option -dir $terraindir/sugarloaf sugarloaf.in 

$dem2fds $option -overlap -dir $terraindir/tower tower.in 

$dem2fds $option -dir $terraindir/trails trails.in 
$dem2fds $option -dir $terraindir/trails trails2.in 
