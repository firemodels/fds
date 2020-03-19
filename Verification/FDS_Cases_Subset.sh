#!/bin/bash

$QFDS -d Complex_Geometry geom_obst.fds

$QFDS -d Aerosols aerosol_scrubbing.fds
$QFDS -d Aerosols aerosol_agglomeration.fds
$QFDS -d Aerosols propane_flame_deposition.fds
$QFDS -d Aerosols soot_oxidation_wall.fds

$QFDS -d Controls activate_vents.fds
$QFDS -d Controls control_test_2.fds
$QFDS -d Controls hrr_freeze.fds

$QFDS -d Detectors beam_detector.fds
$QFDS -d Detectors smoke_detector.fds

$QFDS -d HVAC ashrae7_quadratic.fds
$QFDS -d HVAC HVAC_aircoil.fds
$QFDS -d HVAC HVAC_energy_pressure.fds
$QFDS -d HVAC HVAC_filter.fds
$QFDS -d HVAC HVAC_mass_transport_combine.fds
$QFDS -d HVAC leak_test.fds

$QFDS -d Sprinklers_and_Sprays bucket_test_1.fds
$QFDS -d Sprinklers_and_Sprays cascade.fds
$QFDS -d Sprinklers_and_Sprays terminal_velocity_dt_0_001.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_1.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_6.fds
$QFDS -d Sprinklers_and_Sprays water_fuel_sprays.fds

