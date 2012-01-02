#!/bin/bash

$RUNFDSMPI 5 Flowfields simple_duct
$RUNFDSMPI 8 Pressure_Solver duct_flow
$RUNFDSMPI 5 Pressure_Solver hallways

