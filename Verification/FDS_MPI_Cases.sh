#!/bin/bash

$RUNFDSMPI 5 Flowfields simple_duct

$RUNFDSMPI 4 Pressure_Solver dancing_eddies_tight
$RUNFDSMPI 4 Pressure_Solver dancing_eddies_default
$RUNFDSMPI 8 Pressure_Solver duct_flow
$RUNFDSMPI 5 Pressure_Solver hallways
$RUNFDSMPI 8 Pressure_Solver scarc2d_bicg_8mesh
$RUNFDSMPI 8 Pressure_Solver scarc2d_fft_8mesh
$RUNFDSMPI 8 Pressure_Solver scarc2d_cg_8mesh
$RUNFDSMPI 8 Pressure_Solver scarc2d_gmg_8mesh


