% McDermott
% 5-28-2009
% master_verification_script.m
%
% If you author a section in the verification guide, create a script
% that generates all the graphics (in pdf format) in that section and 
% store the script in the 'scripts' directory.  For example, wall_model.m
% creates all the pdfs for the section on the Werner and Wengle wall
% model.  Also, add your script to the master list below.
%
% The most important script is dataplot.m.  To utilize this script, add the
% appropriate parameters to a 'd' line in verification_data_config_matlab.csv.
%
% If you create your own script (your plots are not generated with dataplot.m),
% please do not forget to add the SVN number to the plot.  Examples for how
% this is done can be found in any of the 'Special cases' scripts below dataplot.

close all
clear all

addpath 'scripts'

% Scripts to run prior to dataplot

radiation_box;            display('radiation_box complete...')
radiation_plane_layer;    display('radiation_plane_layer complete...')
ns2d;                     display('ns2d complete...')
wall_internal_radiation;  display('wall_internal_radiation complete...')
ashrae_7;                 display('ashrae_7 complete...')
flame_species;            display('flame_species complete...')

% The main plotting routine is dataplot

cfil = [pwd,'/verification_data_config_matlab.csv'];
vdir = [pwd,'/../../Verification/'];
plotdir = [pwd,'/../../Manuals/'];

[saved_data,drange] = dataplot(cfil,vdir,plotdir);

% Special cases
 
run scripts/turb_model
run scripts/wall_model
run scripts/pyrolysis
run scripts/birch_tga
run scripts/water_ice_water
run scripts/pcm_slab
run scripts/pulsating
run scripts/compression_wave
run scripts/plate_view_factor
run scripts/low_flux_hot_gas_filling
run scripts/terminal_velocity_convergence
run scripts/flat_fire_comparison
run scripts/fluid_part

display('verification scripts completed successfully!')
