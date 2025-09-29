#!$FIREMODELS/fds/.github/fds_python_env/bin/python

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# # Scripts to run prior to dataplot

# print("catchpole_spread_rates...");       subprocess.run(["python","./scripts/catchpole_spread_rates.py"])
# print("NIST_deposition_gauge...");        subprocess.run(["python","./scripts/NIST_deposition_gauge.py"])
# print("flame_height...");                 subprocess.run(["python","./scripts/flame_height.py"])
# print("NIST_RSE...");                     subprocess.run(["python","./scripts/NIST_RSE.py"])
# print("sippola_aerosol_deposition...");   subprocess.run(["python","./scripts/sippola_aerosol_deposition.py"])
# print("layer_height...");                 subprocess.run(["python","./scripts/layer_height.py"])
# print("NIST_NRC_Corner_Effects...");      subprocess.run(["python","./scripts/NIST_NRC_Corner_Effects.py"])
# # print("fm_data_center...");               subprocess.run(["python","./scripts/fm_data_center.py"])
# print("LNG_Dispersion...");               subprocess.run(["python","./scripts/LNG_Dispersion.py"])
# print("LNG_wind_profiles...");            subprocess.run(["python","./scripts/LNG_wind_profiles.py"])
# print("FM_Vertical_Wall_Flames...");      subprocess.run(["python","./scripts/FM_Vertical_Wall_Flames.py"])
# print("umd_line_burner_process...");      subprocess.run(["python","./scripts/umd_line_burner_process.py"])
# # print("Askervein_Hill...");               subprocess.run(["python","./scripts/Askervein_Hill.py"])
# print("UWO_Wind_Tunnel...");              subprocess.run(["python","./scripts/UWO_Wind_Tunnel.py"])
# print("FM_Burner...");                    subprocess.run(["python","./scripts/FM_Burner.py"])
# print("Crown_Fires...");                  subprocess.run(["python","./scripts/Crown_Fires.py"])
# print("Ranz_Marshall...");                subprocess.run(["python","./scripts/Ranz_Marshall.py"])
# print("Phoenix_LNG_Fires...");            subprocess.run(["python","./scripts/Phoenix_LNG_Fires.py"])

# Dataplot and scatplot options

# Statistics output options

# # Run dataplot and scatplot scripts

# fdsplotlib.dataplot(config_filename='../Matlab/FDS_validation_dataplot_inputs.csv',
#                     expdir='../../../exp/',
#                     cmpdir='../../../out/',
#                     pltdir='../../Manuals/',
#                     close_figs=True,
#                     verbose=True,
#                     plot_range=["all"]) # see notes below on plot_range

# Special cases

print("FM_FPRF_Datacenter...");    subprocess.run(["python","./scripts/FM_FPRF_Datacenter.py"])
print("McCaffrey_Plume...");       subprocess.run(["python","./scripts/McCaffrey_Plume.py"])
print("Sandia_Pools...");          subprocess.run(["python","./scripts/Sandia_Pools.py"])
print("TUS_Facade...");            subprocess.run(["python","./scripts/TUS_Facade_contours.py"])
print("USFS_Deep_Fuel_Beds...");   subprocess.run(["python","./scripts/USFS_Deep_Fuel_Beds.py"])
print("Wu_Bakar_Tunnels...");      subprocess.run(["python","./scripts/Wu_Bakar_Tunnels.py"])

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

