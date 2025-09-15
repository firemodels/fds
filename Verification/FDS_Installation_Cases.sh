#!/bin/bash

$QFDS -d Controls activate_vents.fds
$QFDS -d Energy_Budget energy_budget_cold_walls.fds
$QFDS -d Energy_Budget energy_budget_tmix.fds
$QFDS -d Fires spray_burner.fds
$QFDS -d Heat_Transfer heat_conduction_kc.fds
$QFDS -d HVAC ashrae7_table.fds
$QFDS -d Pressure_Effects isentropic2.fds
$QFDS -p 8 -d Pressure_Solver duct_flow.fds
$QFDS -p 8 -d Pressure_Solver duct_flow_uglmat.fds
$QFDS -d Pyrolysis enthalpy.fds
$QFDS -d Pyrolysis shrink_swell.fds
$QFDS -d Radiation radiation_shield.fds
$QFDS -d Species reactionrate_arrhenius_2order_1step.fds
$QFDS -d Species reactionrate_EDC_1step_CH4_nonmix.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U50_N16.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_1.fds

