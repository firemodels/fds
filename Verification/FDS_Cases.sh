#!/bin/bash

$QFDS -p 4 -d Adaptive_Mesh_Refinement random_meshes.fds
# $QFDS -p 2 -d Adaptive_Mesh_Refinement ns2d_16_emb_1to1_refinement.fds
# $QFDS -p 2 -d Adaptive_Mesh_Refinement ns2d_16_emb_1to2_refinement.fds
# $QFDS -p 13 -d Adaptive_Mesh_Refinement ns2d_16_int_1to2_refinement.fds

$QFDS -d Aerosols aerosol_agglomeration.fds
$QFDS -d Aerosols aerosol_agglomeration_2.fds
$QFDS -d Aerosols aerosol_gravitational_deposition.fds
$QFDS -d Aerosols aerosol_gravitational_deposition_2.fds
$QFDS -d Aerosols aerosol_thermophoretic_deposition.fds
$QFDS -d Aerosols aerosol_thermophoretic_deposition_2.fds
$QFDS -d Aerosols aerosol_turbulent_deposition.fds
$QFDS -d Aerosols aerosol_scrubbing.fds
$QFDS -d Aerosols propane_flame_deposition.fds
$QFDS -d Aerosols propane_flame_deposition_gravitational.fds
$QFDS -d Aerosols propane_flame_deposition_none.fds
$QFDS -d Aerosols propane_flame_deposition_thermophoretic.fds
$QFDS -d Aerosols propane_flame_deposition_turbulent.fds
$QFDS -d Aerosols soot_oxidation_wall.fds

$QFDS -p 36 -d Atmospheric_Effects atmospheric_boundary_layer_1.fds
$QFDS -p 36 -d Atmospheric_Effects atmospheric_boundary_layer_2.fds
$QFDS -p 36 -d Atmospheric_Effects atmospheric_boundary_layer_3.fds
$QFDS -p 36 -d Atmospheric_Effects atmospheric_boundary_layer_4.fds
$QFDS -d Atmospheric_Effects lee_waves.fds
$QFDS -d Atmospheric_Effects stack_effect.fds
$QFDS -d Atmospheric_Effects lapse_rate.fds
$QFDS -p 8 -d Atmospheric_Effects wind_example_5.fds
$QFDS -p 8 -d Atmospheric_Effects wind_example_10.fds
$QFDS -p 32 -d Atmospheric_Effects wind_example_32.fds
$QFDS -d Atmospheric_Effects MO_velocity_profile_stable.fds
$QFDS -d Atmospheric_Effects MO_velocity_profile_unstable.fds

$QFDS -d Chemistry ign_delay_Methane_grimech30_T900K_Phi0p6.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T900K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T900K_Phi1p4.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1000K_Phi0p6.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1000K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1000K_Phi1p4.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1100K_Phi0p6.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1100K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1100K_Phi1p4.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1200K_Phi0p6.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1200K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Methane_grimech30_T1200K_Phi1p4.fds
$QFDS -d Chemistry ign_delay_Methane_TianfengLu_T1100K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Methane_Smooke_T1100K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Ethylene_TianfengLu_T1100K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Propane_USC_T1100K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_Propane_Z66_T1100K_Phi1p0.fds
$QFDS -d Chemistry ign_delay_nHeptane_Chalmers_T1100K_Phi1p0.fds

$QFDS -p 24 -d Chemistry load_bal_methane_smooke_serial.fds
$QFDS -p 24 -d Chemistry load_bal_methane_smooke_parallel.fds
$QFDS -p 6  -d Chemistry load_bal_propane_2step_arrhenius_serial.fds
$QFDS -p 6  -d Chemistry load_bal_propane_2step_arrhenius_parallel.fds
$QFDS -p 6  -d Chemistry load_bal_methane_2step_fast_serial.fds
$QFDS -p 6  -d Chemistry load_bal_methane_2step_fast_parallel.fds

$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_Zeta1p0.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_Zeta0p75.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_Zeta0p5.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_Zeta0p25.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_Zeta0p0.fds
$QFDS -d Chemistry EDC_MultiCFDStep_Methane_grimech30_Zeta1p0.fds
$QFDS -d Chemistry EDC_MultiCFDStep_Methane_grimech30_Zeta0p75.fds
$QFDS -d Chemistry EDC_MultiCFDStep_Methane_grimech30_Zeta0p5.fds
$QFDS -d Chemistry EDC_MultiCFDStep_Methane_grimech30_Zeta0p25.fds
$QFDS -d Chemistry EDC_MultiCFDStep_Methane_grimech30_Zeta0p0.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_taumix0p1.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_taumix0p01.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_taumix0p001.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_taumix0p0001.fds
$QFDS -d Chemistry EDC_OneCFDStep_Methane_grimech30_taumix0p00001.fds

$QFDS -p 24 -d Chemistry EDC_load_bal_methane_smooke_serial.fds
$QFDS -p 24 -d Chemistry EDC_load_bal_methane_smooke_parallel.fds

$QFDS -d Complex_Geometry geom_channel.fds
$QFDS -d Complex_Geometry geom_channel2.fds
$QFDS -d Complex_Geometry geom_channel_tmp.fds
$QFDS -d Complex_Geometry geom_channel_tmp2.fds
$QFDS -d Complex_Geometry geom_extruded_poly.fds
$QFDS -d Complex_Geometry geom_hrrpuv_init.fds
$QFDS -d Complex_Geometry geom_intersection.fds
$QFDS -d Complex_Geometry geom_bad_inverted_normals.fds
$QFDS -p 4 -d Complex_Geometry geom_particle_cascade_2.fds
$QFDS -d Complex_Geometry geom_simple.fds
$QFDS -d Complex_Geometry geom_obst.fds
$QFDS -d Complex_Geometry geom_part_init.fds
$QFDS -d Complex_Geometry geom_sphere1a.fds
$QFDS -d Complex_Geometry geom_sphere1c.fds
$QFDS -d Complex_Geometry geom_sphere1e.fds
$QFDS -d Complex_Geometry geom_sphere3b.fds
$QFDS -d Complex_Geometry geom_sphere3d.fds
$QFDS -d Complex_Geometry geom_sphere3f.fds
$QFDS -d Complex_Geometry geom_sphere1b.fds
$QFDS -d Complex_Geometry geom_sphere1d.fds
$QFDS -d Complex_Geometry geom_sphere1f.fds
$QFDS -d Complex_Geometry geom_sphere3a.fds
$QFDS -d Complex_Geometry geom_sphere3c.fds
$QFDS -d Complex_Geometry geom_sphere3e.fds
$QFDS -d Complex_Geometry geom_terrain.fds
$QFDS -d Complex_Geometry geom_terrain2.fds
$QFDS -d Complex_Geometry geom_texture.fds
$QFDS -d Complex_Geometry geom_texture2.fds
$QFDS -d Complex_Geometry geom_texture3a.fds
$QFDS -d Complex_Geometry geom_texture3b.fds
$QFDS -d Complex_Geometry geom_texture4a.fds
$QFDS -d Complex_Geometry geom_texture4b.fds
$QFDS -d Complex_Geometry sphere_helium_1mesh.fds
$QFDS -p 3 -d Complex_Geometry sphere_helium_3meshes.fds
$QFDS -p 3 -d Complex_Geometry sphere_helium_conserve_3meshes.fds
$QFDS -d Complex_Geometry sphere_radiate.fds
$QFDS -d Complex_Geometry sphere_leak.fds
$QFDS -d Complex_Geometry saad_CC_explicit_512_cfl_p25.fds
$QFDS -d Complex_Geometry saad_CC_explicit_512_cfl_p125.fds
$QFDS -d Complex_Geometry saad_CC_explicit_512_cfl_p0625.fds
$QFDS -d Complex_Geometry shunn3_32_cc_exp_chm.fds
$QFDS -d Complex_Geometry shunn3_64_cc_exp_chm.fds
$QFDS -d Complex_Geometry shunn3_128_cc_exp_chm.fds
$QFDS -d Complex_Geometry shunn3_256_cc_exp_chm.fds
$QFDS -o 2 -d Complex_Geometry shunn3_320_cc_exp_chm.fds
$QFDS -d Complex_Geometry shunn3_32_cc_exp_gdv.fds
$QFDS -d Complex_Geometry shunn3_64_cc_exp_gdv.fds
$QFDS -d Complex_Geometry shunn3_128_cc_exp_gdv.fds
$QFDS -d Complex_Geometry shunn3_256_cc_exp_gdv.fds
$QFDS -o 2 -d Complex_Geometry shunn3_320_cc_exp_gdv.fds
$QFDS -d Complex_Geometry t34_scaling.fds

$QFDS -d Complex_Geometry geom_poiseuille_N10a_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N20a_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N40a_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N80a_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N10na_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N20na_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N40na_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N80na_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N10nah_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N20nah_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N40nah_theta0_stm.fds
$QFDS -d Complex_Geometry geom_poiseuille_N80nah_theta0_stm.fds

$QFDS -d Complex_Geometry geom_bad_inconsistent_normals.fds
$QFDS -d Complex_Geometry geom_bad_non_manifold_edge.fds
$QFDS -d Complex_Geometry geom_bad_non_manifold_vert.fds
$QFDS -d Complex_Geometry geom_bad_open_surface.fds
$QFDS -d Complex_Geometry geom_self_intersection.fds

$QFDS -d Complex_Geometry rotated_cube_0deg_32_stm.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_32_obs.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_64_stm.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_64_obs.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_128_stm.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_128_obs.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_256_stm.fds
$QFDS -d Complex_Geometry rotated_cube_0deg_256_obs.fds
$QFDS -o 2 -d Complex_Geometry rotated_cube_0deg_320_stm.fds
$QFDS -o 2 -d Complex_Geometry rotated_cube_0deg_320_obs.fds
$QFDS -d Complex_Geometry rotated_cube_27deg_32_stm.fds
$QFDS -d Complex_Geometry rotated_cube_27deg_64_stm.fds
$QFDS -d Complex_Geometry rotated_cube_27deg_128_stm.fds
$QFDS -d Complex_Geometry rotated_cube_27deg_256_stm.fds
$QFDS -o 2 -d Complex_Geometry rotated_cube_27deg_320_stm.fds
$QFDS -d Complex_Geometry rotated_cube_45deg_32_stm.fds
$QFDS -d Complex_Geometry rotated_cube_45deg_64_stm.fds
$QFDS -d Complex_Geometry rotated_cube_45deg_128_stm.fds
$QFDS -d Complex_Geometry rotated_cube_45deg_256_stm.fds
$QFDS -o 2 -d Complex_Geometry rotated_cube_45deg_320_stm.fds

#$QFDS -d Complex_Geometry zero_thick_roof.fds
$QFDS -d Complex_Geometry cone_1mesh.fds
$QFDS -d Complex_Geometry geom_mass_file_test.fds
$QFDS -p 5 -d Complex_Geometry geom_stretched_grid.fds
$QFDS -p 3 -d Complex_Geometry thin_object_mass.fds

$QFDS -d Controls activate_vents.fds
$QFDS -d Controls control_test.fds
$QFDS -d Controls control_test_2.fds
$QFDS -d Controls create_remove.fds
$QFDS -d Controls cycle_test.fds
$QFDS -d Controls device_test.fds
$QFDS -d Controls hrr_freeze.fds
$QFDS -d Controls rms_cov_corr.fds
$QFDS -d Controls rms_example.fds
$QFDS -d Controls special_time_ramps.fds
$QFDS -p 4 -d Controls spreading_fire.fds
# if [ "$QUEUE" == "" ]; then
#     QUEUE='batch'
# fi
# sbatch --partition=$QUEUE $FIREMODELS/fds/Verification/Controls/ext_heartbeat_std_curve_slurm.sh
# $QFDS -d Controls ext_heartbeat_std_curve.fds

$QFDS -d Detectors aspiration_detector.fds
$QFDS -p 8 -d Detectors beam_detector.fds
$QFDS -d Detectors objects_dynamic.fds
$QFDS -d Detectors objects_static.fds
$QFDS -d Detectors smoke_detector.fds

$QFDS -d Energy_Budget energy_budget_adiabatic_walls.fds
$QFDS -d Energy_Budget energy_budget_adiabatic_two_fuels.fds
$QFDS -d Energy_Budget energy_budget_cold_walls.fds
$QFDS -d Energy_Budget energy_budget_combustion.fds
$QFDS -d Energy_Budget energy_budget_dns_100.fds
$QFDS -d Energy_Budget energy_budget_particles.fds
$QFDS -d Energy_Budget energy_budget_solid.fds
$QFDS -d Energy_Budget energy_budget_tmix.fds
$QFDS -d Energy_Budget test_hrr_2d_cyl.fds
$QFDS -d Energy_Budget test_hrr_3d.fds

$QFDS -d Extinction extinction_1.fds
$QFDS -d Extinction extinction_2.fds

$QFDS -p 4 -d Fires box_burn_away1.fds
$QFDS -d Fires box_burn_away2.fds
$QFDS -d Fires box_burn_away3.fds
$QFDS -d Fires box_burn_away4.fds
$QFDS -d Fires box_burn_away5.fds
$QFDS -d Fires box_burn_away6.fds
$QFDS -d Fires box_burn_away7.fds
$QFDS -d Fires box_burn_away8.fds
$QFDS -p 4 -d Fires box_burn_away9.fds
$QFDS -p 4 -d Fires box_burn_away10.fds
$QFDS -p 4 -d Fires box_burn_away11.fds
$QFDS -d Fires box_burn_away_2D.fds
$QFDS -d Fires box_burn_away_2D_residue.fds
$QFDS -p 4 -d Fires couch.fds
$QFDS -d Fires fire_const_gamma.fds
$QFDS -d Fires spray_burner.fds
$QFDS -d Fires HoC_Ideal.fds
$QFDS -d Fires HoC_NonIdeal.fds
$QFDS -d Fires simple_test.fds
$QFDS -d Fires tmp_lower_limit_simple.fds
$QFDS -d Fires tmp_lower_limit_default.fds
$QFDS -d Fires tmp_lower_limit_dt_p001.fds
$QFDS -d Fires tunnel_const_gamma.fds
$QFDS -d Fires tunnel_linear_cp.fds
$QFDS -p 8 -d Fires circular_burner.fds

$QFDS -d Flowfields divergence_test_1.fds
$QFDS -d Flowfields divergence_test_2.fds
$QFDS -p 4 -d Flowfields divergence_test_3.fds
$QFDS -d Flowfields cyl_test_1.fds
$QFDS -d Flowfields cyl_test_2.fds
$QFDS -d Flowfields cyl_test_3.fds
$QFDS -d Flowfields cyl_test_4.fds
$QFDS -d Flowfields gas_filling.fds
$QFDS -d Flowfields helium_1d_isothermal.fds
$QFDS -d Flowfields helium_1d_const_gamma.fds
$QFDS -d Flowfields helium_2d_isothermal.fds
$QFDS -d Flowfields helium_air_jet_floor.fds
$QFDS -d Flowfields hole.fds
$QFDS -d Flowfields no_hole.fds
$QFDS -d Flowfields jet_fan.fds
$QFDS -d Flowfields symmetry_test.fds
$QFDS -d Flowfields symmetry_test_2.fds
$QFDS -d Flowfields tangential_velocity.fds
$QFDS -d Flowfields velocity_bc_test.fds
$QFDS -d Flowfields blasius_16.fds
$QFDS -d Flowfields blasius_32.fds
$QFDS -d Flowfields blasius_64.fds
$QFDS -d Flowfields mwtest_cfl.fds
$QFDS -d Flowfields mass_heat_wall_device_test.fds
$QFDS -d Flowfields mass_heat_wall_device_test_2.fds
$QFDS -d Flowfields particle_offgas_1.fds
$QFDS -d Flowfields particle_offgas_2.fds
$QFDS -d Flowfields particle_offgas_3.fds
$QFDS -d Flowfields particle_offgas_4.fds
$QFDS -d Flowfields species_conservation_1.fds
$QFDS -d Flowfields species_conservation_2.fds
$QFDS -d Flowfields species_conservation_3.fds
$QFDS -d Flowfields species_conservation_4.fds
$QFDS -d Flowfields hot_layer_360.fds
$QFDS -d Flowfields realizable_mass_fractions.fds
$QFDS -p 4 -d Flowfields parabolic_profile.fds
$QFDS -p 5 -d Flowfields simple_duct.fds
$QFDS -p 8 -d Flowfields symmetry_test_mpi.fds
$QFDS -p 8 -d Flowfields volume_flow_1.fds
$QFDS -p 8 -d Flowfields volume_flow_2.fds

$QFDS -d Heat_Transfer adiabatic_con_flux.fds
$QFDS -d Heat_Transfer adiabatic_net_flux.fds
$QFDS -p 4 -d Heat_Transfer checkerboard.fds
$QFDS -d Heat_Transfer convective_cooling.fds
$QFDS -d Heat_Transfer convective_cooling_p1.fds
$QFDS -d Heat_Transfer convective_cooling_p05.fds
$QFDS -d Heat_Transfer convective_cooling_p025.fds
$QFDS -d Heat_Transfer convective_cooling_p01.fds
$QFDS -d Heat_Transfer convective_cooling_p005.fds
$QFDS -d Heat_Transfer convective_cooling_p0025.fds
$QFDS -d Heat_Transfer convective_cooling_p00125.fds
$QFDS -d Heat_Transfer heat_conduction_a.fds
$QFDS -d Heat_Transfer heat_conduction_b.fds
$QFDS -d Heat_Transfer heat_conduction_c.fds
$QFDS -d Heat_Transfer heat_conduction_d.fds
$QFDS -d Heat_Transfer heat_conduction_kc.fds
$QFDS -d Heat_Transfer insulated_steel_column.fds
$QFDS -d Heat_Transfer insulated_steel_pipe.fds
$QFDS -d Heat_Transfer insulated_steel_pipe_2d.fds
$QFDS -d Heat_Transfer insulated_steel_plate.fds
$QFDS -d Heat_Transfer ht1d_pile.fds
$QFDS -d Heat_Transfer ht3d_pile.fds
$QFDS -d Heat_Transfer ht3d_beam_heating_1.fds
$QFDS -d Heat_Transfer ht3d_beam_heating_2.fds
$QFDS -p 8 -d Heat_Transfer ht3d_demo.fds
$QFDS -d Heat_Transfer ht3d_energy_conservation.fds
$QFDS -d Heat_Transfer ht3d_energy_conservation_2.fds
$QFDS -d Heat_Transfer ht3d_energy_conservation_3.fds
$QFDS -p 8 -d Heat_Transfer ht3d_energy_conservation_4.fds
$QFDS -p 8 -d Heat_Transfer ht3d_energy_conservation_5.fds
$QFDS -p 8 -d Heat_Transfer ht3d_energy_conservation_6.fds
$QFDS -p 8 -d Heat_Transfer ht3d_energy_conservation_7.fds
$QFDS -p 36 -d Heat_Transfer ht3d_energy_conservation_8.fds
$QFDS -p 8 -d Heat_Transfer ht3d_energy_conservation_9.fds
$QFDS -d Heat_Transfer ht3d_ibeam.fds
$QFDS -d Heat_Transfer ht3d_mass_conservation.fds
$QFDS -d Heat_Transfer ht3d_mass_conservation_2.fds
$QFDS -p 2 -d Heat_Transfer ht3d_network.fds
$QFDS -d Heat_Transfer ht3d_slab.fds
$QFDS -d Heat_Transfer ht3d_sphere_24.fds
$QFDS -p 4 -d Heat_Transfer ht3d_sphere_48.fds
$QFDS -p 64 -d Heat_Transfer ht3d_sphere_96.fds
$QFDS -p 4 -d Heat_Transfer back_wall_test.fds
$QFDS -p 3 -d Heat_Transfer back_wall_test_2.fds
$QFDS -d Heat_Transfer SFPE_Case_1.fds
$QFDS -d Heat_Transfer SFPE_Case_2.fds
$QFDS -d Heat_Transfer SFPE_Case_6.fds
$QFDS -d Heat_Transfer SFPE_Case_7.fds
$QFDS -d Heat_Transfer SFPE_Case_8.fds
$QFDS -d Heat_Transfer SFPE_Case_9.fds
$QFDS -d Heat_Transfer SFPE_Case_13.fds
$QFDS -p 8 -d Heat_Transfer SFPE_Case_16.fds
$QFDS -d Heat_Transfer thermocouple_time_constant.fds

$QFDS -d HVAC ashrae7_fixed_flow.fds
$QFDS -d HVAC ashrae7_quadratic.fds
$QFDS -d HVAC ashrae7_table.fds
$QFDS -d HVAC door_crack.fds
$QFDS -d HVAC fan_test.fds
$QFDS -d HVAC qfan_test.fds
$QFDS -d HVAC qfan_multi.fds
$QFDS -d HVAC HVAC_aircoil.fds
$QFDS -d HVAC HVAC_damper.fds
$QFDS -d HVAC HVAC_energy_pressure.fds
$QFDS -d HVAC HVAC_filter.fds
$QFDS -d HVAC HVAC_flow_loss.fds
$QFDS -d HVAC HVAC_mass_conservation.fds
$QFDS -d HVAC HVAC_mass_transport.fds
$QFDS -d HVAC HVAC_mass_transport_1.fds
$QFDS -d HVAC HVAC_mass_transport_2.fds
$QFDS -d HVAC HVAC_mass_transport_3.fds
$QFDS -d HVAC HVAC_mass_transport_4.fds
$QFDS -d HVAC HVAC_mass_transport_branch.fds
$QFDS -d HVAC HVAC_mass_transport_combine.fds
$QFDS -d HVAC HVAC_mass_transport_conservation_1.fds
$QFDS -d HVAC HVAC_mass_transport_conservation_2.fds
$QFDS -d HVAC HVAC_mass_transport_conservation_3.fds
$QFDS -d HVAC HVAC_mass_transport_conv_0020.fds
$QFDS -d HVAC HVAC_mass_transport_conv_0040.fds
$QFDS -d HVAC HVAC_mass_transport_conv_0080.fds
$QFDS -d HVAC HVAC_mass_transport_conv_0160.fds
$QFDS -d HVAC HVAC_mass_transport_conv_0320.fds
$QFDS -d HVAC HVAC_mass_transport_energy.fds
$QFDS -d HVAC HVAC_tee_loss_1.fds
$QFDS -d HVAC HVAC_tee_loss_2.fds
$QFDS -d HVAC leak_test.fds
$QFDS -d HVAC leak_test_2.fds
$QFDS -d HVAC leak_enthalpy.fds
$QFDS -p 2 -d HVAC leak_test_3.fds
$QFDS -d HVAC leak_test_4.fds
$QFDS -d HVAC HVAC_leak_exponent.fds
$QFDS -p 8 -d HVAC HVAC_leak_exponent_2.fds

$QFDS -d Miscellaneous layer_1mesh.fds
$QFDS -p 3 -d Miscellaneous layer_4mesh.fds
$QFDS -d Miscellaneous pyramid.fds
$QFDS -d Miscellaneous mesh_transformation.fds
$QFDS -d Miscellaneous obst_sphere.fds
$QFDS -d Miscellaneous obst_cylinder.fds
$QFDS -d Miscellaneous obst_cone.fds
$QFDS -d Miscellaneous obst_multi.fds
$QFDS -d Miscellaneous obst_single.fds
$QFDS -d Miscellaneous obst_rotbox.fds
$QFDS -d Miscellaneous obst_sphere_mass_flux.fds
$QFDS -d Miscellaneous obst_cylinder_mass_flux.fds
$QFDS -d Miscellaneous obst_cone_mass_flux.fds
$QFDS -d Miscellaneous obst_box_mass_flux.fds
$QFDS -d Miscellaneous part_orientation_ramp.fds
$QFDS -p 2 -d Miscellaneous part_path_ramp_jog.fds
$QFDS -d Miscellaneous devc_interpolation_velocity.fds
$QFDS -d Miscellaneous devc_interpolation_temperature.fds
$QFDS -d Miscellaneous init_overlap.fds
$QFDS -p 2 -d Miscellaneous external_test.fds
$QFDS -d Miscellaneous ramp_time.fds

$QFDS -d NS_Analytical_Solution ns2d_16.fds
$QFDS -d NS_Analytical_Solution ns2d_16_nupt1.fds
$QFDS -d NS_Analytical_Solution ns2d_32.fds
$QFDS -d NS_Analytical_Solution ns2d_32_nupt1.fds
$QFDS -d NS_Analytical_Solution ns2d_64.fds
$QFDS -d NS_Analytical_Solution ns2d_64_nupt1.fds
$QFDS -d NS_Analytical_Solution ns2d_8.fds
$QFDS -d NS_Analytical_Solution ns2d_8_nupt1.fds
$QFDS -d NS_Analytical_Solution vort2d_40.fds
$QFDS -d NS_Analytical_Solution vort2d_80.fds
$QFDS -d NS_Analytical_Solution vort2d_160.fds
$QFDS -d NS_Analytical_Solution vort2d_320.fds

$QFDS -d Pressure_Effects isentropic.fds
$QFDS -d Pressure_Effects isentropic2.fds
$QFDS -d Pressure_Effects pressure_boundary.fds
$QFDS -d Pressure_Effects pressure_rise.fds
$QFDS -d Pressure_Effects zone_break_fast.fds
$QFDS -d Pressure_Effects zone_break_slow.fds
$QFDS -d Pressure_Effects zone_break_fast_ulmat.fds
$QFDS -d Pressure_Effects zone_break_slow_ulmat.fds
$QFDS -p 4 -d Pressure_Effects zone_break_fast_uglmat.fds
$QFDS -p 2 -d Pressure_Effects zone_shape.fds
$QFDS -p 8 -d Pressure_Effects zone_shape_2.fds
$QFDS -d Pressure_Effects obst_coarse_fine_interface.fds

$QFDS -d Pressure_Solver dancing_eddies_1mesh.fds
$QFDS -d Pressure_Solver opening_ulmat.fds
$QFDS -p 5 -d Pressure_Solver dancing_eddies_embed.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_ulmat.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_ulmat_hypre.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_uglmat.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_uglmat_hypre.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_tight.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_tight_no_precon.fds
$QFDS -p 4 -d Pressure_Solver dancing_eddies_default.fds
$QFDS -p 8 -d Pressure_Solver duct_flow.fds
$QFDS -p 8 -d Pressure_Solver duct_flow_ulmat.fds
$QFDS -p 8 -d Pressure_Solver duct_flow_ulmat_hypre.fds
$QFDS -p 8 -d Pressure_Solver duct_flow_uglmat.fds
$QFDS -p 8 -d Pressure_Solver duct_flow_uglmat_hypre.fds
$QFDS -p 5 -d Pressure_Solver hallways.fds
$QFDS -p 8 -d Pressure_Solver tunnel_demo.fds
$QFDS -p 8 -d Pressure_Solver tunnel_demo_glmat.fds
$QFDS -p 8 -d Pressure_Solver pressure_iteration2d_default.fds
$QFDS -p 8 -d Pressure_Solver pressure_iteration2d_uglmat.fds
$QFDS -p 8 -d Pressure_Solver pressure_iteration3d_default.fds
$QFDS -p 8 -d Pressure_Solver pressure_iteration3d_uglmat.fds
$QFDS -p 16 -d Pressure_Solver random_obstructions_fft.fds
$QFDS -p 4 -d Pressure_Solver stairwell.fds
$QFDS -p 2 -d Pressure_Solver ulmat_2zone.fds
$QFDS -p 4 -d Pressure_Solver obst_activation_default.fds
$QFDS -p 4 -d Pressure_Solver obst_activation_ulmat.fds
$QFDS -p 4 -d Pressure_Solver obst_activation_default_gases.fds

$QFDS -d Pyrolysis anca-couce-fig1_10K.fds
$QFDS -d Pyrolysis anca-couce-fig1_5K.fds
$QFDS -d Pyrolysis anca-couce-fig1_2p5K.fds
$QFDS -d Pyrolysis anca-couce-fig2_10K.fds
$QFDS -d Pyrolysis anca-couce-fig2_5K.fds
$QFDS -d Pyrolysis anca-couce-fig2_2p5K.fds
$QFDS -d Pyrolysis cable_11_insulation_mcc.fds
$QFDS -d Pyrolysis cable_23_insulation_mcc.fds
$QFDS -d Pyrolysis cable_701_insulation_mcc.fds
$QFDS -d Pyrolysis cable_11_jacket_mcc.fds
$QFDS -d Pyrolysis cable_23_jacket_mcc.fds
$QFDS -d Pyrolysis cable_701_jacket_mcc.fds
$QFDS -d Pyrolysis cell_burn_away.fds
$QFDS -d Pyrolysis cellulose_TGA_1KPM_Air.fds
$QFDS -d Pyrolysis cellulose_TGA_3KPM_Air.fds
$QFDS -d Pyrolysis cellulose_TGA_5KPM_Air.fds
$QFDS -d Pyrolysis birch_tga_1step_2.fds
$QFDS -d Pyrolysis birch_tga_1step_20.fds
$QFDS -d Pyrolysis enthalpy.fds
$QFDS -d Pyrolysis ice_cube.fds
$QFDS -d Pyrolysis liquid_mixture.fds
$QFDS -p 8 -d Pyrolysis methanol_evaporation.fds
$QFDS -d Pyrolysis part_baking_soda_420K.fds
$QFDS -d Pyrolysis part_baking_soda_450K.fds
$QFDS -d Pyrolysis part_baking_soda_500K.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp03_1C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp03_3C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp13_1C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp13_3C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp14_1C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp14_3C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp15_1C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp15_3C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp16_1C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp16_3C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp17_1C.fds
$QFDS -d Pyrolysis pine_wood_TGA_exp17_3C.fds
$QFDS -d Pyrolysis pyrolysis_1.fds
$QFDS -d Pyrolysis pyrolysis_2.fds
$QFDS -d Pyrolysis specified_hrr.fds
$QFDS -d Pyrolysis shrink_swell.fds
$QFDS -d Pyrolysis spyro_cone_demo.fds
$QFDS -d Pyrolysis spyro_cone_demo_2.fds
$QFDS -d Pyrolysis spyro_cone_demo_3.fds
$QFDS -d Pyrolysis surf_mass_vent_liquid_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_liquid_fuel_nonconforming.fds
$QFDS -d Pyrolysis surf_mass_part_char_cart_fuel.fds
$QFDS -d Pyrolysis surf_mass_part_char_cart_gas.fds
$QFDS -d Pyrolysis surf_mass_part_char_cyl_fuel.fds
$QFDS -d Pyrolysis surf_mass_part_char_cyl_gas.fds
$QFDS -d Pyrolysis surf_mass_part_char_spher_fuel.fds
$QFDS -d Pyrolysis surf_mass_part_char_spher_gas.fds
$QFDS -d Pyrolysis surf_mass_part_nonchar_cart_fuel.fds
$QFDS -d Pyrolysis surf_mass_part_nonchar_cart_gas.fds
$QFDS -d Pyrolysis surf_mass_part_nonchar_cyl_fuel.fds
$QFDS -d Pyrolysis surf_mass_part_nonchar_cyl_gas.fds
$QFDS -d Pyrolysis surf_mass_part_nonchar_spher_fuel.fds
$QFDS -d Pyrolysis surf_mass_part_nonchar_spher_gas.fds
$QFDS -d Pyrolysis surf_mass_vent_char_cart_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_char_cart_gas.fds
$QFDS -d Pyrolysis surf_mass_vent_char_cyl_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_char_cyl_gas.fds
$QFDS -d Pyrolysis surf_mass_vent_char_spher_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_char_spher_gas.fds
$QFDS -d Pyrolysis surf_mass_vent_nonchar_cart_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_nonchar_cart_gas.fds
$QFDS -d Pyrolysis surf_mass_vent_nonchar_cyl_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_nonchar_cyl_gas.fds
$QFDS -d Pyrolysis surf_mass_vent_nonchar_spher_fuel.fds
$QFDS -d Pyrolysis surf_mass_vent_nonchar_spher_gas.fds
$QFDS -d Pyrolysis surf_mass_part_specified.fds
$QFDS -d Pyrolysis surf_mass_two_species_cart.fds
$QFDS -d Pyrolysis surf_mass_two_species_cyl.fds
$QFDS -d Pyrolysis surf_mass_two_species_spher.fds
$QFDS -p 2 -d Pyrolysis tga_analysis.fds
$QFDS -d Pyrolysis tga_sample.fds
$QFDS -d Pyrolysis two_step_solid_reaction.fds
$QFDS -d Pyrolysis water_pool.fds
$QFDS -d Pyrolysis matl_e_cons_1.fds
$QFDS -d Pyrolysis matl_e_cons_2.fds
$QFDS -d Pyrolysis matl_e_cons_3.fds
$QFDS -d Pyrolysis matl_e_cons_4.fds
$QFDS -d Pyrolysis matl_e_cons_5.fds
$QFDS -d Pyrolysis matl_e_cons_6.fds
$QFDS -d Pyrolysis matl_e_cons_7.fds
$QFDS -d Pyrolysis matl_e_cons_8.fds
$QFDS -d Pyrolysis matl_e_cons_9.fds

$QFDS -d Radiation adiabatic_surface_temperature.fds
$QFDS -d Radiation droplet_absorption_cart.fds
$QFDS -d Radiation droplet_absorption_cyl.fds
$QFDS -d Radiation emissivity.fds
$QFDS -p 2 -d Radiation geom_rad.fds
$QFDS -p 8 -d Radiation geom_rad_2.fds
$QFDS -d Radiation hot_spheres.fds
$QFDS -d Radiation part_attenuation.fds
$QFDS -d Radiation plate_view_factor_2D_30.fds
$QFDS -d Radiation plate_view_factor_2D_60.fds
$QFDS -d Radiation plate_view_factor_2D_100.fds
$QFDS -d Radiation plate_view_factor_cart_30.fds
$QFDS -d Radiation plate_view_factor_cart_60.fds
$QFDS -d Radiation plate_view_factor_cart_100.fds
$QFDS -d Radiation plate_view_factor_cyl_30.fds
$QFDS -d Radiation plate_view_factor_cyl_60.fds
$QFDS -d Radiation plate_view_factor_cyl_100.fds
$QFDS -d Radiation plate_view_factor_ibm_30.fds
$QFDS -d Radiation plate_view_factor_ibm_60.fds
$QFDS -d Radiation plate_view_factor_ibm_100.fds
$QFDS -d Radiation radiating_polygon_square_20.fds
$QFDS -d Radiation radiating_polygon_square_40.fds
$QFDS -d Radiation radiating_polygon_square_80.fds
$QFDS -d Radiation radiation_box_100_1000.fds
$QFDS -d Radiation radiation_box_100__100.fds
$QFDS -d Radiation radiation_box_100_2000.fds
$QFDS -d Radiation radiation_box_100__300.fds
$QFDS -d Radiation radiation_box_100___50.fds
$QFDS -d Radiation radiation_box__20_1000.fds
$QFDS -d Radiation radiation_box__20__100.fds
$QFDS -d Radiation radiation_box__20_2000.fds
$QFDS -d Radiation radiation_box__20__300.fds
$QFDS -d Radiation radiation_box__20___50.fds
$QFDS -p 3 -d Radiation radiation_gas_panel.fds
$QFDS -d Radiation radiation_plane_layer_1_1.fds
$QFDS -d Radiation radiation_plane_layer_1_2.fds
$QFDS -d Radiation radiation_plane_layer_1_3.fds
$QFDS -d Radiation radiation_plane_layer_1_4.fds
$QFDS -d Radiation radiation_plane_layer_1_5.fds
$QFDS -d Radiation radiation_plane_layer_2_1.fds
$QFDS -d Radiation radiation_plane_layer_2_2.fds
$QFDS -d Radiation radiation_plane_layer_2_3.fds
$QFDS -d Radiation radiation_plane_layer_2_4.fds
$QFDS -d Radiation radiation_plane_layer_2_5.fds
$QFDS -d Radiation radiation_plane_layer_3_1.fds
$QFDS -d Radiation radiation_plane_layer_3_2.fds
$QFDS -d Radiation radiation_plane_layer_3_3.fds
$QFDS -d Radiation radiation_plane_layer_3_4.fds
$QFDS -d Radiation radiation_plane_layer_3_5.fds
$QFDS -d Radiation radiation_plane_layer_4_1.fds
$QFDS -d Radiation radiation_plane_layer_4_2.fds
$QFDS -d Radiation radiation_plane_layer_4_3.fds
$QFDS -d Radiation radiation_plane_layer_4_4.fds
$QFDS -d Radiation radiation_plane_layer_4_5.fds
$QFDS -d Radiation radiation_plane_layer_5_1.fds
$QFDS -d Radiation radiation_plane_layer_5_2.fds
$QFDS -d Radiation radiation_plane_layer_5_3.fds
$QFDS -d Radiation radiation_plane_layer_5_4.fds
$QFDS -d Radiation radiation_plane_layer_5_5.fds
$QFDS -d Radiation radiation_plane_layer_6_1.fds
$QFDS -d Radiation radiation_plane_layer_6_2.fds
$QFDS -d Radiation radiation_plane_layer_6_3.fds
$QFDS -d Radiation radiation_plane_layer_6_4.fds
$QFDS -d Radiation radiation_plane_layer_6_5.fds
$QFDS -d Radiation radiation_shield.fds
$QFDS -d Radiation target_test.fds
$QFDS -d Radiation thermocouples.fds
$QFDS -d Radiation TC_heating.fds
$QFDS -d Radiation TC_view_factor.fds
$QFDS -d Radiation wall_internal_radiation.fds
$QFDS -d Radiation net_rad_1000_cgeom.fds
$QFDS -d Radiation net_rad_1000_obst.fds

$QFDS -p 4 -d Restart device_restart_a.fds
$QFDS -p 4 -d Restart device_restart_base_case.fds
$QFDS -d Restart restart_test1a.fds
$QFDS -d Restart geom_restart_a.fds
$QFDS -d Restart geom_restart_base_case.fds
$QFDS -d Restart geom_ls_restart_a.fds
$QFDS -d Restart geom_ls_restart_base_case.fds
$QFDS -d Restart clocks_restart_a.fds
$QFDS -d Restart clocks_restart_base_case.fds
$QFDS -p 4 -d Restart restart_ulmat_a.fds
$QFDS -p 8 -d Restart csvf_restart_a.fds

$QFDS -d Species burke_schumann.fds
$QFDS -d Species FED_FIC.fds
$QFDS -d Species FED_FIC_SMIX.fds
$QFDS -d Species FED_CO_HCN.fds
$QFDS -d Species FED_moving.fds
$QFDS -d Species mass_flux_comparison.fds
$QFDS -d Species methane_flame_simple.fds
$QFDS -d Species methane_flame_simple_2.fds
$QFDS -d Species methane_flame_primitive.fds
$QFDS -d Species methane_flame_primitive_2.fds
$QFDS -d Species methane_flame_lumped.fds
$QFDS -d Species methane_flame_lumped_fuel.fds
$QFDS -d Species methane_flame_lumped_ox.fds
$QFDS -d Species propane_flame_2reac.fds
$QFDS -d Species propane_flame_2reac_simple.fds
$QFDS -d Species multiple_reac_hrrpua.fds
$QFDS -d Species multiple_reac_n_simple.fds
$QFDS -d Species reactionrate_arrhenius_0order_1step.fds
$QFDS -d Species reactionrate_arrhenius_0order_1step_cvode.fds
$QFDS -d Species reactionrate_arrhenius_2order_1step.fds
$QFDS -d Species reactionrate_arrhenius_2order_1step_cvode.fds
$QFDS -d Species reactionrate_arrhenius_1p75order_2step.fds
$QFDS -d Species reactionrate_arrhenius_1p75order_2step_cvode.fds
$QFDS -d Species reactionrate_arrhenius_1p75order_2stepr.fds
$QFDS -d Species reactionrate_arrhenius_1p75order_2stepr_cvode.fds
$QFDS -d Species reactionrate_arrhenius_equilibrium.fds
$QFDS -d Species reactionrate_arrhenius_equilibrium_cvode.fds
$QFDS -d Species reactionrate_arrhenius_jones_lindstedt.fds
$QFDS -d Species reactionrate_arrhenius_jones_lindstedt_cvode.fds
$QFDS -d Species reactionrate_EDC_1step_CH4_nonmix.fds
$QFDS -d Species reactionrate_EDC_flim_1step_C3H8.fds
$QFDS -d Species reactionrate_EDC_flim_1step_CH4.fds
$QFDS -d Species reactionrate_EDC_flim_2step.fds
$QFDS -d Species reactionrate_EDC_O2lim_1step.fds
$QFDS -d Species reactionrate_EDC_O2lim_2fuel_prim.fds
$QFDS -d Species reactionrate_EDC_O2lim_2fuel_lump.fds
$QFDS -d Species reactionrate_fast_slow.fds
$QFDS -d Species reactionrate_lumped_two_air.fds
$QFDS -d Species reactionrate_lumped_two_air_2.fds
$QFDS -d Species reactionrate_series_reaction.fds
$QFDS -d Species pvc_combustion.fds
$QFDS -d Species hrrpuv_reac_simple.fds
$QFDS -d Species hrrpuv_reac_extinction.fds
$QFDS -d Species hrrpuv_reac_single.fds
$QFDS -d Species hrrpuv_reac_parallel.fds
$QFDS -d Species hrrpuv_reac_parallel_2.fds
$QFDS -d Species hrrpuv_reac_series.fds
$QFDS -d Species hrrpuv_reac_soot.fds
$QFDS -d Species hrrpuv_reac_arrhenius.fds
$QFDS -d Species ramp_chi_r.fds
$QFDS -d Species bound_test_1.fds
$QFDS -d Species bound_test_2.fds
$QFDS -d Species species_props.fds
$QFDS -d Species condensation_1.fds
$QFDS -d Species condensation_2.fds
$QFDS -d Species condensation_3.fds
$QFDS -d Species wall_cond.fds
$QFDS -d Species humidity.fds
$QFDS -d Species mass_flux_wall_yindex.fds
$QFDS -d Species mass_flux_wall_zindex.fds
$QFDS -d Species mass_balance_gas_volume.fds
$QFDS -p 8 -d Species mass_balance_reac.fds
$QFDS -d Species mass_balance_reac_2.fds
$QFDS -p 3 -d Species favre_test.fds
$QFDS -p 4 -d Species 1_step_2_step_compare.fds

$QFDS -d Sprinklers_and_Sprays activate_sprinklers.fds
$QFDS -d Sprinklers_and_Sprays adiabatic_surf_temp_spray.fds
$QFDS -p 4 -d Sprinklers_and_Sprays bucket_test_1.fds
$QFDS -d Sprinklers_and_Sprays bucket_test_2.fds
$QFDS -d Sprinklers_and_Sprays bucket_test_3.fds
$QFDS -d Sprinklers_and_Sprays bucket_test_4.fds
$QFDS -d Sprinklers_and_Sprays cannon_ball.fds
$QFDS -d Sprinklers_and_Sprays cascade.fds
$QFDS -p 6 -d Sprinklers_and_Sprays cascadempi.fds
$QFDS -p 2 -d Sprinklers_and_Sprays e_coefficient.fds
$QFDS -d Sprinklers_and_Sprays droplet_distributions.fds
$QFDS -d Sprinklers_and_Sprays droplet_distributions_2.fds
$QFDS -d Sprinklers_and_Sprays flow_rate.fds
$QFDS -d Sprinklers_and_Sprays flow_rate_2.fds
$QFDS -p 2 -d Sprinklers_and_Sprays geom_sprk_mass.fds
$QFDS -d Sprinklers_and_Sprays particle_anisotropic_radi.fds
$QFDS -d Sprinklers_and_Sprays particle_isotropic_radi.fds
$QFDS -d Sprinklers_and_Sprays particle_colors.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U10_N16.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U50_N16.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U100_N16.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U50_N1600.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U100_N1600.fds
$QFDS -d Sprinklers_and_Sprays particle_drag_U150_N1600.fds
$QFDS -d Sprinklers_and_Sprays particle_flux.fds
$QFDS -d Sprinklers_and_Sprays particle_heating_convection.fds
$QFDS -d Sprinklers_and_Sprays particle_heating_radiation.fds
$QFDS -d Sprinklers_and_Sprays sphere_drag_1.fds
$QFDS -d Sprinklers_and_Sprays sphere_drag_2.fds
$QFDS -d Sprinklers_and_Sprays terminal_velocity_dt_1_0.fds
$QFDS -d Sprinklers_and_Sprays terminal_velocity_dt_0_1.fds
$QFDS -d Sprinklers_and_Sprays terminal_velocity_dt_0_01.fds
$QFDS -d Sprinklers_and_Sprays terminal_velocity_dt_0_001.fds
$QFDS -d Sprinklers_and_Sprays terminal_velocity_dt_0_0001.fds
$QFDS -d Sprinklers_and_Sprays flat_fire.fds
$QFDS -d Sprinklers_and_Sprays fluid_part_mom_x.fds
$QFDS -d Sprinklers_and_Sprays fluid_part_mom_y.fds
$QFDS -d Sprinklers_and_Sprays fluid_part_mom_z.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_1.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_2.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_3.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_4.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_4_np100.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_5.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_6.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_6a.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_7.fds
$QFDS -d Sprinklers_and_Sprays water_evaporation_8.fds
$QFDS -d Sprinklers_and_Sprays water_fuel_sprays.fds
$QFDS -d Sprinklers_and_Sprays screen_drag_1.fds
$QFDS -d Sprinklers_and_Sprays screen_drag_2.fds
$QFDS -d Sprinklers_and_Sprays porous_media.fds
$QFDS -d Sprinklers_and_Sprays water_evap_1_const_gamma.fds

$QFDS -d Scalar_Analytical_Solution pulsating_FL0_16.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL0_32.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL0_64.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL0_128.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL2_16.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL2_32.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL2_64.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL2_128.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL4_16.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL4_32.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL4_64.fds
$QFDS -d Scalar_Analytical_Solution pulsating_FL4_128.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL0_16.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL0_32.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL0_64.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL0_128.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL2_16.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL2_32.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL2_64.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL2_128.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL4_16.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL4_32.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL4_64.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL4_128.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL5_16.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL5_32.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL5_64.fds
$QFDS -d Scalar_Analytical_Solution compression_wave_FL5_128.fds
$QFDS -d Scalar_Analytical_Solution move_slug.fds
$QFDS -d Scalar_Analytical_Solution move_slug_fl1.fds
$QFDS -d Scalar_Analytical_Solution shunn3_32.fds
$QFDS -d Scalar_Analytical_Solution shunn3_64.fds
$QFDS -d Scalar_Analytical_Solution shunn3_128.fds
$QFDS -d Scalar_Analytical_Solution shunn3_256.fds
$QFDS -d Scalar_Analytical_Solution shunn3_512.fds
$QFDS -d Scalar_Analytical_Solution shunn3_256_cfl_1.fds
$QFDS -d Scalar_Analytical_Solution shunn3_256_cfl_p5.fds
$QFDS -d Scalar_Analytical_Solution shunn3_256_cfl_p25.fds
$QFDS -d Scalar_Analytical_Solution shunn3_256_cfl_p125.fds
$QFDS -d Scalar_Analytical_Solution shunn3_256_cfl_p0625.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_1.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p5.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p25.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p125.fds
$QFDS -d Scalar_Analytical_Solution saad_512_cfl_p0625.fds
$QFDS -d Scalar_Analytical_Solution shunn3_FavreZ_32.fds
$QFDS -d Scalar_Analytical_Solution shunn3_FavreZ_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_32.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_128.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_256.fds
$QFDS -p 4 -d Scalar_Analytical_Solution shunn3_4mesh_512.fds
$QFDS      -d Scalar_Analytical_Solution soborot_godunov_square_wave_16.fds
$QFDS      -d Scalar_Analytical_Solution soborot_godunov_square_wave_32.fds
$QFDS      -d Scalar_Analytical_Solution soborot_godunov_square_wave_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution soborot_godunov_square_wave_128.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_square_wave_16.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_square_wave_32.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_square_wave_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution soborot_superbee_square_wave_128.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_square_wave_128_1mesh.fds
$QFDS      -d Scalar_Analytical_Solution soborot_charm_square_wave_16.fds
$QFDS      -d Scalar_Analytical_Solution soborot_charm_square_wave_32.fds
$QFDS      -d Scalar_Analytical_Solution soborot_charm_square_wave_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution soborot_charm_square_wave_128.fds
$QFDS      -d Scalar_Analytical_Solution soborot_charm_cos_wave_16.fds
$QFDS      -d Scalar_Analytical_Solution soborot_charm_cos_wave_32.fds
$QFDS      -d Scalar_Analytical_Solution soborot_charm_cos_wave_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution soborot_charm_cos_wave_128.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_cos_wave_16.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_cos_wave_32.fds
$QFDS      -d Scalar_Analytical_Solution soborot_superbee_cos_wave_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution soborot_superbee_cos_wave_128.fds
$QFDS      -d Scalar_Analytical_Solution soborot_mp5_cos_wave_16.fds
$QFDS      -d Scalar_Analytical_Solution soborot_mp5_cos_wave_32.fds
$QFDS      -d Scalar_Analytical_Solution soborot_mp5_cos_wave_64.fds
$QFDS -p 4 -d Scalar_Analytical_Solution soborot_mp5_cos_wave_128.fds

$QFDS -o 1 -d Thread_Check race_test_1.fds
$QFDS -o 4 -d Thread_Check race_test_4.fds

$QFDS -o 1 -d Timing_Benchmarks openmp_test64a.fds
$QFDS -o 2 -d Timing_Benchmarks openmp_test64b.fds
$QFDS -o 3 -d Timing_Benchmarks openmp_test64c.fds
$QFDS -o 4 -d Timing_Benchmarks openmp_test64d.fds
$QFDS -o 5 -d Timing_Benchmarks openmp_test64e.fds
$QFDS -o 6 -d Timing_Benchmarks openmp_test64f.fds
$QFDS -o 7 -d Timing_Benchmarks openmp_test64g.fds
$QFDS -o 8 -d Timing_Benchmarks openmp_test64h.fds
$QFDS -o 1 -d Timing_Benchmarks openmp_test128a.fds
$QFDS -o 2 -d Timing_Benchmarks openmp_test128b.fds
$QFDS -o 3 -d Timing_Benchmarks openmp_test128c.fds
$QFDS -o 4 -d Timing_Benchmarks openmp_test128d.fds
$QFDS -o 5 -d Timing_Benchmarks openmp_test128e.fds
$QFDS -o 6 -d Timing_Benchmarks openmp_test128f.fds
$QFDS -o 7 -d Timing_Benchmarks openmp_test128g.fds
$QFDS -o 8 -d Timing_Benchmarks openmp_test128h.fds

$QFDS -d Turbulence csmag0_32.fds
$QFDS -d Turbulence csmag_32.fds
$QFDS -d Turbulence csmag_64.fds
$QFDS -d Turbulence dsmag_32.fds
$QFDS -d Turbulence dsmag_64.fds
$QFDS -d Turbulence mu0_32.fds
$QFDS -d Turbulence deardorff_32.fds
$QFDS -d Turbulence deardorff_64.fds
$QFDS -d Turbulence vreman_32.fds
$QFDS -d Turbulence vreman_64.fds
$QFDS -d Turbulence wale_32.fds
$QFDS -d Turbulence wale_64.fds
$QFDS -d Turbulence yplus_8.fds
$QFDS -d Turbulence yplus_16.fds
$QFDS -d Turbulence yplus_32.fds
$QFDS -d Turbulence sem_flat_leddy_p2.fds
$QFDS -d Turbulence sem_par_leddy_p2.fds
$QFDS -d Turbulence sem_atm_leddy_p2.fds
$QFDS -d Turbulence sem_ramp_leddy_p2.fds
$QFDS -d Turbulence sem_open_wind.fds
$QFDS -d Turbulence ramp_prof_u_z.fds

$QFDS -d WUI Bova_1a.fds
$QFDS -d WUI Bova_1b.fds
$QFDS -d WUI Bova_4a.fds
$QFDS -d WUI bulk_density_file.fds
$QFDS -d WUI char_oxidation_1.fds
$QFDS -p 2 -d WUI char_oxidation_2.fds
$QFDS -d WUI dragon_5a.fds
$QFDS -p 2 -d WUI ground_vegetation_conv.fds
$QFDS -p 3 -d WUI ground_vegetation_drag.fds
$QFDS -p 2 -d WUI ground_vegetation_load.fds
$QFDS -p 2 -d WUI ground_vegetation_radi.fds
$QFDS -d WUI hot_rods.fds
$QFDS -d WUI level_set_fuel_model_1.fds
$QFDS -d WUI LS_ellipse_0ms_00deg.fds
$QFDS -d WUI LS_ellipse_5ms_00deg.fds
$QFDS -d WUI LS_ellipse_0ms_30deg.fds
$QFDS -d WUI LS_ellipse_5ms_30deg.fds
$QFDS -d WUI LS_wind_ramp_lin.fds
$QFDS -d WUI LS_wind_ramp_quad.fds
$QFDS -d WUI LS4_ember_ignition.fds
$QFDS -d WUI LS4_ember_yield.fds
$QFDS -d WUI Morvan_TGA.fds
$QFDS -d WUI Morvan_TGA_2.fds
$QFDS -d WUI Needle_TGA_A.fds
$QFDS -d WUI Needle_TGA_B.fds
$QFDS -d WUI Needle_TGA_C.fds
$QFDS -d WUI Needle_TGA_D.fds
$QFDS -d WUI Needle_TGA_E.fds
$QFDS -d WUI Needle_TGA_F.fds
$QFDS -d WUI Needle_TGA_G.fds
$QFDS -d WUI Needle_TGA_H.fds
$QFDS -d WUI part_drag_prof_ux.fds
$QFDS -d WUI part_drag_prof_uy.fds
$QFDS -d WUI part_drag_prof_uz.fds
$QFDS -d WUI part_drag_prof_vx.fds
$QFDS -d WUI part_drag_prof_vy.fds
$QFDS -d WUI part_drag_prof_vz.fds
$QFDS -d WUI part_drag_prof_wx.fds
$QFDS -d WUI part_drag_prof_wy.fds
$QFDS -d WUI part_drag_prof_wz.fds
$QFDS -d WUI part_drag_stretched.fds
$QFDS -d WUI pine_needles.fds
$QFDS -d WUI radiation_gas-veg_consistency_gas.fds
$QFDS -d WUI radiation_gas-veg_consistency_veg.fds
$QFDS -d WUI random_walk_1.fds
$QFDS -d WUI random_walk_2.fds
$QFDS -d WUI tree_shapes.fds
$QFDS -p 6 -d WUI vegetation_absorb.fds
$QFDS -d WUI vegetation_drag_1.fds
$QFDS -d WUI vegetation_drag_2.fds
