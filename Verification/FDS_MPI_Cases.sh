#!/bin/bash

$RUNFDSMPI 5 Flowfields simple_duct
$RUNFDSMPI 8 Pressure_Solver duct_flow
$RUNFDSMPI 5 Pressure_Solver hallways
$RUNFDSMPI 16 Pressure_Solver scarc2d_bicg_16mesh
$RUNFDSMPI 16 Pressure_Solver scarc2d_cggmg_16mesh
$RUNFDSMPI 16 Pressure_Solver scarc2d_bicggmg_16mesh

