#!$FIREMODELS/fds/.github/fds_python_env/bin/python

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)

# Scripts to run prior to dataplot

# print("catchpole_spread_rates...");   subprocess.run(["python","./scripts/catchpole_spread_rates.py"])
# print("NIST_deposition_gauge...");   subprocess.run(["python","./scripts/NIST_deposition_gauge.py"])
# print("flame_height...");   subprocess.run(["python","./scripts/flame_height.py"])
# print("NIST_RSE...");   subprocess.run(["python","./scripts/NIST_RSE.py"])
# print("sippola_aerosol_deposition...");   subprocess.run(["python","./scripts/sippola_aerosol_deposition.py"])
# print("layer_height...");   subprocess.run(["python","./scripts/layer_height.py"])
# print("NIST_NRC_Corner_Effects...");   subprocess.run(["python","./scripts/NIST_NRC_Corner_Effects.py"])
# # print("fm_data_center...");   subprocess.run(["python","./scripts/fm_data_center.py"])
# print("LNG_Dispersion...");   subprocess.run(["python","./scripts/LNG_Dispersion.py"])
# print("LNG_wind_profiles...");   subprocess.run(["python","./scripts/LNG_wind_profiles.py"])

# Dataplot and scatplot options

# Statistics output options

# # Run dataplot and scatplot scripts

# fdsplotlib.dataplot(config_filename='../Matlab/FDS_validation_dataplot_inputs.csv',
#                     expdir='../../../exp/',
#                     cmpdir='../../../out/',
#                     pltdir='../../Manuals/',
#                     close_figs=True,
#                     verbose=True,
#                     plot_list=['all'])

# Special cases

print("Sandia_Pools...");          subprocess.run(["python","./scripts/Sandia_Pools.py"])
print("TUS_Facade...");            subprocess.run(["python","./scripts/TUS_Facade_contours.py"])
print("Wu_Bakar_Tunnels...");      subprocess.run(["python","./scripts/Wu_Bakar_Tunnels.py"])
print("USFS_Deep_Fuel_Beds...");   subprocess.run(["python","./scripts/USFS_Deep_Fuel_Beds.py"])

print("Python validation scripts completed successfully!")
