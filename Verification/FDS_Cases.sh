#!/bin/bash

# add -A to any case that you wish to be a part of the benchmark timing suite

$QFDS -d Complex_Geometry saad_CC_explicit_512_cfl_p25.fds
$QFDS -d Complex_Geometry saad_CC_explicit_512_cfl_p125.fds
$QFDS -d Complex_Geometry saad_CC_explicit_512_cfl_p0625.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_1.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p5.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p25.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p125.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p0625.fds
