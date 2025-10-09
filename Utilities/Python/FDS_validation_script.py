#!$FIREMODELS/fds/.github/fds_python_env/bin/python

import fdsplotlib
import importlib
import runpy
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# Scripts to run prior to dataplot

print("catchpole_spread_rates...");       runpy.run_path("./scripts/catchpole_spread_rates.py", run_name="__main__")
print("NIST_deposition_gauge...");        runpy.run_path("./scripts/NIST_deposition_gauge.py", run_name="__main__")
print("flame_height...");                 runpy.run_path("./scripts/flame_height.py", run_name="__main__")
print("NIST_RSE...");                     runpy.run_path("./scripts/NIST_RSE.py", run_name="__main__")
print("sippola_aerosol_deposition...");   runpy.run_path("./scripts/sippola_aerosol_deposition.py", run_name="__main__")
print("layer_height...");                 runpy.run_path("./scripts/layer_height.py", run_name="__main__")
print("NIST_NRC_Corner_Effects...");      runpy.run_path("./scripts/NIST_NRC_Corner_Effects.py", run_name="__main__")
# print("fm_data_center...");               runpy.run_path("./scripts/fm_data_center.py", run_name="__main__")
print("LNG_Dispersion...");               runpy.run_path("./scripts/LNG_Dispersion.py", run_name="__main__")
print("LNG_wind_profiles...");            runpy.run_path("./scripts/LNG_wind_profiles.py", run_name="__main__")
print("FM_Vertical_Wall_Flames...");      runpy.run_path("./scripts/FM_Vertical_Wall_Flames.py", run_name="__main__")
print("umd_line_burner_process...");      runpy.run_path("./scripts/umd_line_burner_process.py", run_name="__main__")
# print("Askervein_Hill...");               runpy.run_path("./scripts/Askervein_Hill.py", run_name="__main__")
print("UWO_Wind_Tunnel...");              runpy.run_path("./scripts/UWO_Wind_Tunnel.py", run_name="__main__")
print("FM_Burner...");                    runpy.run_path("./scripts/FM_Burner.py", run_name="__main__")
print("Crown_Fires...");                  runpy.run_path("./scripts/Crown_Fires.py", run_name="__main__")
print("Ranz_Marshall...");                runpy.run_path("./scripts/Ranz_Marshall.py", run_name="__main__")
print("Phoenix_LNG_Fires...");            runpy.run_path("./scripts/Phoenix_LNG_Fires.py", run_name="__main__")
print("Sandia_Plumes_TKE...");            runpy.run_path("./scripts/Sandia_Plumes_TKE.py", run_name="__main__")

# # Dataplot and scatplot options

# Dataplot_Inputs_File = '../Matlab/FDS_validation_dataplot_inputs.csv'
# EXP_Dir = '../../../exp/'
# OUT_Dir = '../../../out/'
# Manuals_Dir = '../../Manuals/'
# Scatterplot_Inputs_File = '../Matlab/FDS_validation_scatterplot_inputs.csv'

# # Statistics output options

# Stats_Output = 'Validation'
# Scatterplot_Dir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/'

# # Run dataplot and scatplot scripts

# saved_data, drange = fdsplotlib.dataplot(config_filename=Dataplot_Inputs_File,
#                                          expdir=EXP_Dir,
#                                          cmpdir=OUT_Dir,
#                                          pltdir=Manuals_Dir,
#                                          close_figs=True,
#                                          verbose=True,
#                                          plot_range=["Sandia Plumes"],
#                                          ) # see notes below on plot_range

# fdsplotlib.scatplot(saved_data,drange,
# 				    Manuals_Dir=Manuals_Dir,
# 				    Scatterplot_Inputs_File=Scatterplot_Inputs_File,
# 				    Stats_Output=Stats_Output,
# 				    Scatterplot_Dir=Scatterplot_Dir,
# 				    verbose=True,
# 				    )


# Special cases
print("Beyler_Hood...");                  runpy.run_path("./scripts/Beyler_Hood.py", run_name="__main__")
print("BRE_LEMTA_Sprays...");             runpy.run_path("./scripts/BRE_LEMTA_Sprays.py", run_name="__main__")
print("FHWA_Tunnel...");                  runpy.run_path("./scripts/FHWA_Tunnel.py", run_name="__main__")
print("FM_FPRF_Datacenter...");           runpy.run_path("./scripts/FM_FPRF_Datacenter.py", run_name="__main__")
print("Heskestad_Flame_Height_2...");     runpy.run_path("./scripts/Heskestad_Flame_Height_2.py", run_name="__main__")
print("McCaffrey_Plume...");              runpy.run_path("./scripts/McCaffrey_Plume.py", run_name="__main__")
print("Memorial_Tunnel...");              runpy.run_path("./scripts/Memorial_Tunnel.py", run_name="__main__")
print("Memorial_Tunnel_2...");            runpy.run_path("./scripts/Memorial_Tunnel_2.py", run_name="__main__")
print("NIST_NRC_Parallel_Panels...");     runpy.run_path("./scripts/NIST_NRC_Parallel_Panels.py", run_name="__main__")
print("Sandia_Plumes...");                runpy.run_path("./scripts/Sandia_Plumes.py", run_name="__main__")
print("Sandia_Pools...");                 runpy.run_path("./scripts/Sandia_Pools.py", run_name="__main__")
print("Theobald_Hose_Stream...");         runpy.run_path("./scripts/Theobald_Hose_Stream.py", run_name="__main__")
print("TUS_Facade...");                   runpy.run_path("./scripts/TUS_Facade_contours.py", run_name="__main__")
print("USFS_Deep_Fuel_Beds...");          runpy.run_path("./scripts/USFS_Deep_Fuel_Beds.py", run_name="__main__")
print("Wu_Bakar_Tunnels...");             runpy.run_path("./scripts/Wu_Bakar_Tunnels.py", run_name="__main__")

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

