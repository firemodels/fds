#!/bin/bash
option=$1 $2 
if [ "$option" == ""]; then
option=-o
fi

dem2fds=dem2fds
#dem2fds=../intel_linux_64/dem2fds_linux_64

terraindir=~/terrain

echo demtest
$dem2fds $option -n -d $terraindir/demtest demtest1.in 
$dem2fds $option -n -d $terraindir/demtest demtest2.in 

echo blodget
$dem2fds $option -d $terraindir/blodget blodget.in 

echo NIST
$dem2fds $option -n -d $terraindir/nist nist.in 

echo tower
$dem2fds $option -d $terraindir/tower tower.in 

echo sugarloaf
$dem2fds $option -n -d $terraindir/sugarloaf sugarloaf.in 

echo trails
$dem2fds $option -d $terraindir/trails trails.in 
