#!$FIREMODELS/fds/.github/fds_python_env/bin/python
# McDermott
# 2 April 2024

import subprocess
import fdsplotlib
import runpy
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)
print("Using:", fdsplotlib.__file__)
# Scripts to run prior to dataplot

# print("ignition_delay...");   runpy.run_path("./scripts/cantera_ignition_delay.py", run_name="__main__")
# print("reaction_rates...");   runpy.run_path("./scripts/cantera_reaction_rates.py", run_name="__main__")
# print("turbulent_batch_reactor...");   runpy.run_path("./scripts/cantera_turbulent_batch_reactor.py", run_name="__main__")
print("ashrae_7...");                            runpy.run_path("./scripts/ashrae_7.py", run_name="__main__")
print('burke_schumann...');                      runpy.run_path("./scripts/burke_schumann.py", run_name="__main__")
print('cat_propane_depo...');                    runpy.run_path("./scripts/cat_propane_depo.py", run_name="__main__")
print('convective_cooling_convergence...');      runpy.run_path("./scripts/convective_cooling_convergence.py", run_name="__main__")
print('flame_species...');                       runpy.run_path("./scripts/flame_species.py", run_name="__main__")
print("ns2d...");                                runpy.run_path("./scripts/ns2d.py", run_name="__main__")
print("radiation_box...");                       runpy.run_path("./scripts/radiation_box.py", run_name="__main__")
print("radiation_plane_layer...");               runpy.run_path("./scripts/radiation_plane_layer.py", run_name="__main__")
print("random_walk_soln...");                    runpy.run_path("./scripts/random_walk_soln.py", run_name="__main__")
print('rms_cov_corr...');                        runpy.run_path("./scripts/rms_cov_corr.py", run_name="__main__")
print('vegetation_absorb...');                   runpy.run_path("./scripts/vegetation_absorb.py", run_name="__main__")
print('water_evap_1_const_gamma...');            runpy.run_path("./scripts/water_evap_1_const_gamma.py", run_name="__main__")

# Dataplot and scatplot options

# Statistics output options

# Run dataplot and scatplot scripts

fdsplotlib.dataplot(config_filename='../Matlab/FDS_verification_dataplot_inputs.csv',
                    expdir='../../Verification/',
                    cmpdir='../../Verification/',
                    pltdir='../../Manuals/',
                    close_figs=True,
                    verbose=True,
                    plot_range=[2,2]) # plot_range[start, end], optionally instead use plot_list['Dataname']

# Special cases

print("atmospheric_boundary_layer...");     runpy.run_path("./scripts/atmospheric_boundary_layer.py", run_name="__main__")
print("blasius...");                        runpy.run_path("./scripts/blasius.py", run_name="__main__")
print('compression_wave...');               runpy.run_path("./scripts/compression_wave.py", run_name="__main__")
print('extinction...');                     runpy.run_path("./scripts/extinction.py", run_name="__main__")
print("fan_curve...");                      runpy.run_path("./scripts/fan_curve.py", run_name="__main__")
print("favre_test...");                     runpy.run_path("./scripts/favre_test.py", run_name="__main__")
print("fds_moody_chart...");                runpy.run_path("./scripts/fds_moody_chart.py", run_name="__main__")
print("fluid_part...");                     runpy.run_path("./scripts/fluid_part.py", run_name="__main__")
print("freecon_sphere...");                 runpy.run_path("./scripts/freecon_sphere.py", run_name="__main__")
print("geom_channel_test...");              runpy.run_path("./scripts/geom_channel_test.py", run_name="__main__")
print("geom_positive_errors...");           runpy.run_path("./scripts/geom_positive_errors.py", run_name="__main__")
print("heated_channel...");                 runpy.run_path("./scripts/heated_channel.py", run_name="__main__")
print('hot_layer_collapse...');             runpy.run_path("./scripts/hot_layer_collapse.py", run_name="__main__")
print("ht3d_sphere...");                    runpy.run_path("./scripts/ht3d_sphere.py", run_name="__main__")
print('hvac_mass_transport...');            runpy.run_path("./scripts/hvac_mass_transport.py", run_name="__main__")
print("jet_decay...");                      runpy.run_path("./scripts/jet_decay.py", run_name="__main__")
print("law_of_the_wall...");                runpy.run_path("./scripts/law_of_the_wall.py", run_name="__main__")
print("level_set_ellipse...");              runpy.run_path("./scripts/level_set_ellipse.py", run_name="__main__")
print("mesh_transformation...");            runpy.run_path("./scripts/mesh_transformation.py", run_name="__main__")
print("mass_balance...");                   runpy.run_path("./scripts/mass_balance.py", run_name="__main__")
print("mass_balance_gas_volume...");        runpy.run_path("./scripts/mass_balance_gas_volume.py", run_name="__main__")
print("mass_balance_reac...");              runpy.run_path("./scripts/mass_balance_reac.py", run_name="__main__")
print("openmp_timing_benchmarks...");       runpy.run_path("./scripts/openmp_timing_benchmarks.py", run_name="__main__")
print("particle_size_distribution...");     runpy.run_path("./scripts/particle_size_distribution.py", run_name="__main__")
print("plate_view_factor...");              runpy.run_path("./scripts/plate_view_factor.py", run_name="__main__")
print("pulsating...");                      runpy.run_path("./scripts/pulsating.py", run_name="__main__")
print("pyrolysis...");                      runpy.run_path("./scripts/pyrolysis.py", run_name="__main__")
print('radiating_polygon...');              runpy.run_path("./scripts/radiating_polygon.py", run_name="__main__")
print("ribbed_channel...");                 runpy.run_path("./scripts/ribbed_channel.py", run_name="__main__")
print('rotcube_cc_mms_error...');           runpy.run_path("./scripts/rotcube_cc_mms_error.py", run_name="__main__")
print("saad_mms...");                       runpy.run_path("./scripts/saad_mms_temporal_error.py", run_name="__main__")
print("scaling_tests...");                  runpy.run_path("./scripts/scaling_tests.py", run_name="__main__")
print("shunn_mms...");                      runpy.run_path("./scripts/shunn_mms.py", run_name="__main__")
print("soborot_mass_transport...");         runpy.run_path("./scripts/soborot_mass_transport.py", run_name="__main__")
print("synthetic_eddy_method...");          runpy.run_path("./scripts/synthetic_eddy_method.py", run_name="__main__")
print("terminal_velocity_convergence...");  runpy.run_path("./scripts/terminal_velocity_convergence.py", run_name="__main__")
print("tree_shapes...");                    runpy.run_path("./scripts/tree_shapes.py", run_name="__main__")
print("turb_model...");                     runpy.run_path("./scripts/turb_model.py", run_name="__main__")
print("vort2d...");                         runpy.run_path("./scripts/vort2d.py", run_name="__main__")
print("wall_internal_radiation...");        runpy.run_path("./scripts/wall_internal_radiation.py", run_name="__main__")
print("yplus...");                          runpy.run_path("./scripts/yplus.py", run_name="__main__")
print('impinging_jet...');                  runpy.run_path("./scripts/impinging_jet.py", run_name="__main__")

print("verification scripts completed successfully!")
