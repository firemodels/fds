#!$FIREMODELS/fds/.github/fds_python_env/bin/python
# McDermott
# 2 April 2024

import subprocess
import fdsplotlib
import matplotlib.pyplot as plt
import runpy
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)
print("Using:", fdsplotlib.__file__)

# If there is an error in one of the sub-scripts, print the message but do not stop the main script.

def safe_run(script_path):
    try:
        runpy.run_path(script_path, run_name="__main__")
        plt.clf()         # Clear the current figure (if any)
        plt.close('all')  # Close all open figure windows
    except Exception as exc:
        print(f"Error in {script_path}: {exc}")

# Scripts to run prior to dataplot

print("ashrae_7...");                            safe_run("./scripts/ashrae_7.py")
print('burke_schumann...');                      safe_run("./scripts/burke_schumann.py")
print('cat_propane_depo...');                    safe_run("./scripts/cat_propane_depo.py")
print('convective_cooling_convergence...');      safe_run("./scripts/convective_cooling_convergence.py")
print('flame_species...');                       safe_run("./scripts/flame_species.py")
print("ns2d...");                                safe_run("./scripts/ns2d.py")
print("radiation_box...");                       safe_run("./scripts/radiation_box.py")
print("radiation_plane_layer...");               safe_run("./scripts/radiation_plane_layer.py")
print("random_walk_soln...");                    safe_run("./scripts/random_walk_soln.py")
print('rms_cov_corr...');                        safe_run("./scripts/rms_cov_corr.py")
print('vegetation_absorb...');                   safe_run("./scripts/vegetation_absorb.py")
print('water_evap_1_const_gamma...');            safe_run("./scripts/water_evap_1_const_gamma.py")

# Dataplot and scatplot options

Dataplot_Inputs_File = 'FDS_verification_dataplot_inputs.csv';
Working_Dir = '../../Verification/';
Manuals_Dir = '../../Manuals/';
Scatterplot_Inputs_File = 'FDS_verification_scatterplot_inputs.csv';

# Statistics output options

Stats_Output = 'Verification'
Scatterplot_Dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/'

# Run dataplot and scatplot scripts

saved_data, drange = fdsplotlib.dataplot(config_filename=Dataplot_Inputs_File,
                                         expdir=Working_Dir,
                                         cmpdir=Working_Dir,
                                         pltdir=Manuals_Dir,
                                         close_figs=True,
                                         verbose=True,
                                         plot_range=["all"]) # see notes below on plot_range

# ----- write saved_data, drange to disk -----
import pickle

# Save
with open("saved_data_verification.pkl", "wb") as f:
    pickle.dump((saved_data, drange), f)

# Later...
with open("saved_data_verification.pkl", "rb") as f:
    saved_data, drange = pickle.load(f)
#---------------------------------------------

fdsplotlib.scatplot(saved_data,drange,
                    Manuals_Dir=Manuals_Dir,
                    Scatterplot_Inputs_File=Scatterplot_Inputs_File,
                    Stats_Output=Stats_Output,
                    Scatterplot_Dir=Scatterplot_Dir,
                    verbose=True,
                    )

# Special cases

print("atmospheric_boundary_layer...");     safe_run("./scripts/atmospheric_boundary_layer.py")
print("blasius...");                        safe_run("./scripts/blasius.py")
print('compression_wave...');               safe_run("./scripts/compression_wave.py")
print('extinction...');                     safe_run("./scripts/extinction.py")
print("fan_curve...");                      safe_run("./scripts/fan_curve.py")
print("favre_test...");                     safe_run("./scripts/favre_test.py")
print("fds_moody_chart...");                safe_run("./scripts/fds_moody_chart.py")
print("fds_timing_stats...");               safe_run("./scripts/fds_timing_stats.py")
print("fluid_part...");                     safe_run("./scripts/fluid_part.py")
print("freecon_sphere...");                 safe_run("./scripts/freecon_sphere.py")
print("geom_channel_test...");              safe_run("./scripts/geom_channel_test.py")
print("geom_positive_errors...");           safe_run("./scripts/geom_positive_errors.py")
print("heated_channel...");                 safe_run("./scripts/heated_channel.py")
print('hot_layer_collapse...');             safe_run("./scripts/hot_layer_collapse.py")
print("ht3d_sphere...");                    safe_run("./scripts/ht3d_sphere.py")
print('hvac_mass_transport...');            safe_run("./scripts/hvac_mass_transport.py")
print("jet_decay...");                      safe_run("./scripts/jet_decay.py")
print("law_of_the_wall...");                safe_run("./scripts/law_of_the_wall.py")
print("level_set_ellipse...");              safe_run("./scripts/level_set_ellipse.py")
print("mesh_transformation...");            safe_run("./scripts/mesh_transformation.py")
print("make_smv_images...");                safe_run("./scripts/make_smv_images.py")
print("mass_balance...");                   safe_run("./scripts/mass_balance.py")
print("mass_balance_gas_volume...");        safe_run("./scripts/mass_balance_gas_volume.py")
print("mass_balance_reac...");              safe_run("./scripts/mass_balance_reac.py")
print("natconh...");                        safe_run("./scripts/natconh.py")
print("natconv...");                        safe_run("./scripts/natconv.py")
print("openmp_timing_benchmarks...");       safe_run("./scripts/openmp_timing_benchmarks.py")
print("particle_size_distribution...");     safe_run("./scripts/particle_size_distribution.py")
print("plate_view_factor...");              safe_run("./scripts/plate_view_factor.py")
print("pulsating...");                      safe_run("./scripts/pulsating.py")
print("pyrolysis...");                      safe_run("./scripts/pyrolysis.py")
print('radiating_polygon...');              safe_run("./scripts/radiating_polygon.py")
print("ribbed_channel...");                 safe_run("./scripts/ribbed_channel.py")
print('rotcube_cc_mms_error...');           safe_run("./scripts/rotcube_cc_mms_error.py")
print("saad_mms...");                       safe_run("./scripts/saad_mms_temporal_error.py")
print("scaling_tests...");                  safe_run("./scripts/scaling_tests.py")
print("shunn_mms...");                      safe_run("./scripts/shunn_mms.py")
print("soborot_mass_transport...");         safe_run("./scripts/soborot_mass_transport.py")
print("synthetic_eddy_method...");          safe_run("./scripts/synthetic_eddy_method.py")
print("terminal_velocity_convergence...");  safe_run("./scripts/terminal_velocity_convergence.py")
print("tree_shapes...");                    safe_run("./scripts/tree_shapes.py")
print("turb_model...");                     safe_run("./scripts/turb_model.py")
print("vort2d...");                         safe_run("./scripts/vort2d.py")
print("wall_internal_radiation...");        safe_run("./scripts/wall_internal_radiation.py")
print("yplus...");                          safe_run("./scripts/yplus.py")
print('impinging_jet...');                  safe_run("./scripts/impinging_jet.py")
print('htc_forced...');                     safe_run("./scripts/htc_forced.py")
print('nat_conv_hot_plate...');             safe_run("./scripts/nat_conv_hot_plate.py")


print("verification scripts completed successfully!")

# ------------------------------
# plot_range usage examples
#
# plot_range lets you select which rows of the config file to process.
# You can mix row numbers, ranges, and Dataname strings:
#
#  1. Single row by number (Spreadsheet-style, including header rows):
#       plot_range = [1995]
#
#  2. Inclusive ranges by "start:stop":
#       plot_range = ["5:9"]        # rows 5 through 9
#
#  3. Open-ended ranges:
#       plot_range = ["1995:"]      # from row 1995 to the end
#
#  4. Named selection by Dataname (case-insensitive):
#       plot_range = ["CSTB Tunnel", "Steckler Compartment"]
#
#  5. Mixed selection:
#       plot_range = [1, 2, "5:9", "CSTB Tunnel", "7000:"]
#
#  6. All rows:
#       plot_range = ["all"]
#
# Notes:
# - Row numbers are 1-based (like Spreadsheet).
# - Ranges are inclusive, e.g. "5:9" means 5,6,7,8,9.
# - "start:" runs to the last row.
# - Strings that are not ranges or "all" are matched to the Dataname column.
# ------------------------------
