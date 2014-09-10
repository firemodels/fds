#!/bin/bash

$QFDS -p 5 -d Flowfields simple_duct.fds

$QFDS -p 4 -d Pressure_Solver dancing_eddies_tight.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_default.fds
$QFDS -p 8 -d Pressure_Solver duct_flow.fds
$QFDS -p 5 -d Pressure_Solver hallways.fds

#$QFDS -p 4 -d WRF wrf_time_ramp.fds
#$QFDS -p 4 -d WRF wrf_prof_ramp.fds
#$QFDS -p 4 -d WRF wrf_time_prof_ramp.fds
