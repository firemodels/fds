#!/bin/bash

$RUNFDSMPI 5 Flowfields simple_duct
$RUNFDSMPI 8 Pressure_Solver duct_flow
$RUNFDSMPI 5 Pressure_Solver hallways
#$RUNFDSMPI 64 Pressure_Solver scarc2d_bicg_64mesh
#$RUNFDSMPI 64 Pressure_Solver scarc2d_cggmg_64mesh
#$RUNFDSMPI 64 Pressure_Solver scarc2d_bicggmg_64mesh

