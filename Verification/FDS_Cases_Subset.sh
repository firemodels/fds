#!/bin/bash

$QFDS -d Complex_Geometry geom_obst.fds
$QFDS -d Complex_Geometry geom_intersection.fds
$QFDS -d Complex_Geometry geom_bad_inverted_normals.fds
$QFDS -d Complex_Geometry geom_simple.fds
$QFDS -d Complex_Geometry sphere_helium_1mesh.fds
$QFDS -d Complex_Geometry sphere_radiate.fds
$QFDS -d Complex_Geometry geom_self_intersection.fds
$QFDS -d Complex_Geometry cone_1mesh.fds

$QFDS -d Aerosols aerosol_scrubbing.fds
$QFDS -d Aerosols aerosol_agglomeration.fds
$QFDS -d Aerosols propane_flame_deposition.fds
$QFDS -d Aerosols soot_oxidation_wall.fds

$QFDS -d Atmospheric_Effects lapse_rate.fds

$QFDS -d Controls activate_vents.fds
$QFDS -d Controls control_test_2.fds
$QFDS -d Controls hrr_freeze.fds

$QFDS -d Detectors beam_detector.fds
$QFDS -d Detectors smoke_detector.fds

$QFDS -d Energy_Budget energy_budget_adiabatic_two_fuels.fds
$QFDS -d Energy_Budget energy_budget_combustion.fds
$QFDS -d Energy_Budget energy_budget_tmix.fds

$QFDS -d Extinction extinction_1.fds
$QFDS -d Extinction extinction_2.fds

$QFDS -d Fires box_burn_away1.fds
$QFDS -d Fires tmp_lower_limit_simple.fds

$QFDS -d Flowfields divergence_test_1.fds
$QFDS -d Flowfields helium_1d_const_gamma.fds
$QFDS -d Flowfields symmetry_test.fds
$QFDS -d Flowfields mwtest_cfl.fds
$QFDS -d Flowfields particle_offgas_1.fds

$QFDS -d Heat_Transfer heat_conduction_a.fds
$QFDS -d Heat_Transfer heat_conduction_b.fds
$QFDS -d Heat_Transfer heat_conduction_c.fds
$QFDS -d Heat_Transfer heat_conduction_d.fds

$QFDS -d HVAC ashrae7_quadratic.fds
$QFDS -d HVAC HVAC_aircoil.fds
$QFDS -d HVAC HVAC_energy_pressure.fds
$QFDS -d HVAC HVAC_filter.fds
$QFDS -d HVAC HVAC_mass_transport_combine.fds
$QFDS -d HVAC leak_test.fds

$QFDS -d NS_Analytical_Solution ns2d_8.fds
$QFDS -d NS_Analytical_Solution ns2d_16.fds
$QFDS -d NS_Analytical_Solution ns2d_32.fds
$QFDS -d NS_Analytical_Solution ns2d_64.fds

$QFDS -d Pyrolysis pyrolysis_1.fds
$QFDS -d Pyrolysis pyrolysis_2.fds
$QFDS -d Pyrolysis two_step_solid_reaction.fds

$QFDS -d Radiation plate_view_factor_cart_30.fds
$QFDS -d Radiation plate_view_factor_cart_60.fds
$QFDS -d Radiation plate_view_factor_cart_100.fds

$QFDS -d Species propane_flame_2reac_simple.fds
$QFDS -d Species hrrpuv_reac_simple.fds

$QFDS -d Sprinklers_and_Sprays bucket_test_1.fds
$QFDS -d Sprinklers_and_Sprays cascade.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_1.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_6.fds
$QFDS -d Sprinklers_and_Sprays water_fuel_sprays.fds

$QFDS -d Turbulence deardorff_32.fds
$QFDS -d Turbulence deardorff_64.fds

$QFDS -d WUI Bova_1a.fds
$QFDS -d WUI Bova_1b.fds
$QFDS -d WUI pine_needles.fds
$QFDS -d WUI vegetation_drag_1.fds
$QFDS -d WUI vegetation_drag_2.fds

