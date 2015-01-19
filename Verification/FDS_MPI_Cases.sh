#!/bin/bash

$QFDS -p 4 -d Adaptive_Mesh_Refinement random_meshes.fds

$QFDS -p 8 -d Fires circular_burner.fds

$QFDS -p 5 -d Flowfields simple_duct.fds
$QFDS -p 5 -d Flowfields simple_duct_2.fds
$QFDS -p 8 -d Flowfields symmetry_test_mpi.fds

$QFDS -p 4 -d Heat_Transfer back_wall_test.fds

$QFDS -p 4 -d Pressure_Solver dancing_eddies_tight.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_tight_overlap.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_default.fds
$QFDS -p 8 -d Pressure_Solver duct_flow.fds
$QFDS -p 5 -d Pressure_Solver hallways.fds
$QFDS -p 8 -d Pressure_Solver scarc2d_bicg_8mesh.fds
$QFDS -p 8 -d Pressure_Solver scarc2d_fft_8mesh.fds
$QFDS -p 8 -d Pressure_Solver scarc2d_cg_8mesh.fds
#$QFDS -p 8 -d Pressure_Solver scarc2d_gmg_8mesh.fds

$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_32.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_128.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_256.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_512.fds

$QFDS -p 4 -d WRF wrf_time_ramp.fds
$QFDS -p 4 -d WRF wrf_prof_ramp.fds
$QFDS -p 4 -d WRF wrf_time_prof_ramp.fds
