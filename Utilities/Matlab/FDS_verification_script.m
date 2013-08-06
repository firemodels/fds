% McDermott
% 5-28-2009
% FDS_verification_script.m
%
% If you author a section in the verification guide, create a script
% that generates all the graphics (in pdf format) in that section and 
% store the script in the 'scripts' directory.  For example, wall_model.m
% creates all the pdfs for the section on the Werner and Wengle wall
% model.  Also, add your script to the master list below.
%
% The most important script is dataplot.m.  To utilize this script, add the
% appropriate parameters to a 'd' line in FDS_verification_dataplot_inputs.csv.
% For debug purposes, you can set the switch_id (first column) of the input
% file to "o" to process only that line.
%
% If you create your own script (your plots are not generated with dataplot.m),
% please do not forget to add the SVN number to the plot.  Examples for how
% this is done can be found in any of the 'Special cases' scripts below dataplot.

close all
clear all

addpath 'scripts'

% Scripts to run prior to dataplot

radiation_box
radiation_plane_layer
ns2d
vort2d
wall_internal_radiation
ashrae_7
flame_species
EDC_species
cutcell_area
cat_propane_depo
burke_schumann

% The main plotting routine is dataplot

cfil = [pwd,'/FDS_verification_dataplot_inputs.csv'];
vdir = [pwd,'/../../Verification/'];
plotdir = [pwd,'/../../Manuals/'];
qfil = [pwd,'/FDS_verification_scatterplot_inputs.csv'];
output_file = [pwd,'/FDS_verification_scatterplot_output.csv'];
stats_output = 1;
statistics_tex_output = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/verification_statistics.tex';

[saved_data,drange] = dataplot(cfil,vdir,plotdir);
scatplot(saved_data,drange,qfil,plotdir,output_file)

% Special cases
 
turb_model
jet_decay
wall_model
pyrolysis
birch_tga
water_ice_water
pcm_slab
pulsating
compression_wave
plate_view_factor
terminal_velocity_convergence
flat_fire_comparison
fluid_part
extinction
fan_curve
mesh_transformation
synthetic_eddy_method

display('verification scripts completed successfully!')
