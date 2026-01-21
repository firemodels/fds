#!$FIREMODELS/fds/.github/fds_python_env/bin/python

import subprocess
import fdsplotlib
import matplotlib.pyplot as plt
import importlib
import runpy
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

print("NIST_deposition_gauge...");        safe_run("./scripts/NIST_deposition_gauge.py")
print("flame_height...");                 safe_run("./scripts/flame_height.py")
print("NIST_RSE...");                     safe_run("./scripts/NIST_RSE.py")
print("sippola_aerosol_deposition...");   safe_run("./scripts/sippola_aerosol_deposition.py")
print("layer_height...");                 safe_run("./scripts/layer_height.py")
print("NIST_NRC_Corner_Effects...");      safe_run("./scripts/NIST_NRC_Corner_Effects.py")
print("LNG_Dispersion...");               safe_run("./scripts/LNG_Dispersion.py")
print("LNG_wind_profiles...");            safe_run("./scripts/LNG_wind_profiles.py")
print("FM_Vertical_Wall_Flames...");      safe_run("./scripts/FM_Vertical_Wall_Flames.py")
print("umd_line_burner_process...");      safe_run("./scripts/umd_line_burner_process.py")
# print("Askervein_Hill...");               safe_run("./scripts/Askervein_Hill.py")
print("UWO_Wind_Tunnel...");              safe_run("./scripts/UWO_Wind_Tunnel.py")
print("FM_Burner...");                    safe_run("./scripts/FM_Burner.py")
print("JIS_Facade_filter...");            safe_run("./scripts/JIS_Facade_filter.py")
print("Crown_Fires...");                  safe_run("./scripts/Crown_Fires.py")
print("Ranz_Marshall...");                safe_run("./scripts/Ranz_Marshall.py")
print("Phoenix_LNG_Fires...");            safe_run("./scripts/Phoenix_LNG_Fires.py")
print("Sandia_Plumes_TKE...");            safe_run("./scripts/Sandia_Plumes_TKE.py")

# Dataplot and scatplot options

Dataplot_Inputs_File = 'FDS_validation_dataplot_inputs.csv'
EXP_Dir = '../../../exp/'
OUT_Dir = '../../../out/'
Manuals_Dir = '../../Manuals/'
Scatterplot_Inputs_File = 'FDS_validation_scatterplot_inputs.csv'

# Statistics output options

Stats_Output = 'Validation'
Scatterplot_Dir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Scatterplots/'

# Run dataplot and scatplot scripts

saved_data, drange = fdsplotlib.dataplot(config_filename=Dataplot_Inputs_File,
                                         expdir=EXP_Dir,
                                         cmpdir=OUT_Dir,
                                         pltdir=Manuals_Dir,
                                         close_figs=True,
                                         verbose=True,
                                         plot_range=["all"],
                                         ) # see notes below on plot_range

# ----- write saved_data, drange to disk -----
import pickle

# Save
with open("saved_data_validation.pkl", "wb") as f:
    pickle.dump((saved_data, drange), f)

# Later...
with open("saved_data_validation.pkl", "rb") as f:
    saved_data, drange = pickle.load(f)
#---------------------------------------------

fdsplotlib.scatplot(saved_data,drange,
                    Manuals_Dir=Manuals_Dir,
                    Scatterplot_Inputs_File=Scatterplot_Inputs_File,
                    Stats_Output=Stats_Output,
                    Scatterplot_Dir=Scatterplot_Dir,
                    verbose=True,
                    )

# Create table of git statistics for FDS Validation Guide

print("validation_git_stats...");         safe_run("./scripts/validation_git_stats.py")

# Special cases

print("Backward_Facing_Step...");         safe_run("./scripts/Backward_Facing_Step.py")
print("Beyler_Hood...");                  safe_run("./scripts/Beyler_Hood.py")
print("BRE_LEMTA_Sprays...");             safe_run("./scripts/BRE_LEMTA_Sprays.py")
print("catchpole_spread_rates...");       safe_run("./scripts/catchpole_spread_rates.py")
print("FHWA_Tunnel...");                  safe_run("./scripts/FHWA_Tunnel.py")
print("FM_FPRF_Datacenter...");           safe_run("./scripts/FM_FPRF_Datacenter.py")
print("Heskestad_Flame_Height_2...");     safe_run("./scripts/Heskestad_Flame_Height_2.py")
print("McCaffrey_Plume...");              safe_run("./scripts/McCaffrey_Plume.py")
print("Memorial_Tunnel...");              safe_run("./scripts/Memorial_Tunnel.py")
print("Memorial_Tunnel_2...");            safe_run("./scripts/Memorial_Tunnel_2.py")
print("NIST_NRC_Parallel_Panels...");     safe_run("./scripts/NIST_NRC_Parallel_Panels.py")
print("NIST_Pool_Fires...");              safe_run("./scripts/NIST_Pool_Fires.py")
print("Sandia_Plumes...");                safe_run("./scripts/Sandia_Plumes.py")
print("Sandia_Pools...");                 safe_run("./scripts/Sandia_Pools.py")
print("Theobald_Hose_Stream...");         safe_run("./scripts/Theobald_Hose_Stream.py")
print("JIS_Facade...");                   safe_run("./scripts/JIS_Facade_contours.py")
print("JIS_Facade_exp...");               safe_run("./scripts/JIS_Facade_exp_contours.py")
print("USFS_Deep_Fuel_Beds...");          safe_run("./scripts/USFS_Deep_Fuel_Beds.py")
print("Wu_Bakar_Tunnels...");             safe_run("./scripts/Wu_Bakar_Tunnels.py")

print("Python validation scripts completed successfully!")

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

