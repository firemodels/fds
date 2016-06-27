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
% please do not forget to add the Git revision number to the plot.  Examples for how
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
cat_propane_depo
burke_schumann
convective_cooling_convergence
random_walk_soln

% Dataplot and scatplot options

Dataplot_Inputs_File = [pwd, '/FDS_verification_dataplot_inputs.csv'];
Working_Dir = [pwd, '/../../FDS/Verification/'];
Manuals_Dir = [pwd, '/../../FDS/Manuals/'];
Scatterplot_Inputs_File = [pwd, '/FDS_verification_scatterplot_inputs.csv'];

% Statistics output options

Stats_Output = 'Verification';
Output_File = [pwd, '/FDS_verification_scatterplot_output.csv'];
Statistics_Tex_Output = [pwd, '/../../FDS/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/verification_statistics.tex'];

% Override the plot style options with NRC 1824 plot options

NRC_Options = false;
Append_To_Scatterplot_Title = '';

% Run dataplot and scatplot scripts

[saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Manuals_Dir);
scatplot(saved_data, drange, ...
         'Manuals_Dir', Manuals_Dir, ...
         'Scatterplot_Inputs_File', Scatterplot_Inputs_File, ...
         'Stats_Output', Stats_Output, ...
         'Output_File', Output_File, ...
         'Statistics_Tex_Output', Statistics_Tex_Output, ...
         'NRC_Options', NRC_Options, ...
         'Append_To_Scatterplot_Title', Append_To_Scatterplot_Title)

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
fluid_part
extinction
fan_curve
mesh_transformation
synthetic_eddy_method
shunn_mms_error
openmp_timing_benchmarks
rms_cov_corr
hot_layer_collapse
radiating_polygon
saad_mms_temporal_error
shunn_mms_temporal_error
scaling_tests

display('verification scripts completed successfully!')
