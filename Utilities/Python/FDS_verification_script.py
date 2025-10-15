#!$FIREMODELS/fds/.github/fds_python_env/bin/python
# McDermott
# 2 April 2024

import subprocess
import fdsplotlib
import importlib
importlib.reload(fdsplotlib) # use for development (while making changes to fdsplotlib.py)
print("Using:", fdsplotlib.__file__)

# Scripts to run prior to dataplot

# print("ignition_delay...");   subprocess.run(["python","./scripts/cantera_ignition_delay.py"])
# print("reaction_rates...");   subprocess.run(["python","./scripts/cantera_reaction_rates.py"])
# print("turbulent_batch_reactor...");   subprocess.run(["python","./scripts/cantera_turbulent_batch_reactor.py"])
print("radiation_box...");   subprocess.run(["python","./scripts/radiation_box.py"])

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

print("blasius...");                        subprocess.run(["python","./scripts/blasius.py"])
print("fan_curve...");                      subprocess.run(["python","./scripts/fan_curve.py"])
print("fds_moody_chart...");                subprocess.run(["python","./scripts/fds_moody_chart.py"])
print("fluid_part...");                     subprocess.run(["python","./scripts/fluid_part.py"])
print("heated_channel...");                 subprocess.run(["python","./scripts/heated_channel.py"])
print("jet_decay...");                      subprocess.run(["python","./scripts/jet_decay.py"])
print("law_of_the_wall...");                subprocess.run(["python","./scripts/law_of_the_wall.py"])
print("level_set_ellipse...");              subprocess.run(["python","./scripts/level_set_ellipse.py"])
print("mesh_transformation...");            subprocess.run(["python","./scripts/mesh_transformation.py"])
print("plate_view_factor...");              subprocess.run(["python","./scripts/plate_view_factor.py"])
print("pulsating...");                      subprocess.run(["python","./scripts/pulsating.py"])
print("pyrolysis...");                      subprocess.run(["python","./scripts/pyrolysis.py"])
print("ribbed_channel...");                 subprocess.run(["python","./scripts/ribbed_channel.py"])
print("saad_mms...");                       subprocess.run(["python","./scripts/saad_mms_temporal_error.py"])
print("shunn_mms...");                      subprocess.run(["python","./scripts/shunn_mms.py"])
print("soborot_mass_transport...");         subprocess.run(["python","./scripts/soborot_mass_transport.py"])
print("synthetic_eddy_method...");          subprocess.run(["python","./scripts/synthetic_eddy_method.py"])
print("terminal_velocity_convergence...");  subprocess.run(["python","./scripts/terminal_velocity_convergence.py"])
print("tree_shapes...");                    subprocess.run(["python","./scripts/tree_shapes.py"])
print("turb_model...");                     subprocess.run(["python","./scripts/turb_model.py"])
print("wall_internal_radiation...");        subprocess.run(["python","./scripts/wall_internal_radiation.py"])
print("yplus...");                          subprocess.run(["python","./scripts/yplus.py"])

print("verification scripts completed successfully!")
