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
water_evap_1_const_gamma

% Dataplot and scatplot options

Dataplot_Inputs_File = 'FDS_verification_dataplot_inputs.csv';
Working_Dir = '../../Verification/';
Manuals_Dir = '../../Manuals/';
Scatterplot_Inputs_File = 'FDS_verification_scatterplot_inputs.csv';

% Statistics output options

Stats_Output = 'Verification';
Output_File = [pwd, '/FDS_verification_scatterplot_output.csv'];
Statistics_Tex_Output = [pwd, '/../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/verification_statistics.tex'];

% Override the plot style options with NRC 1824 plot options

NRC_Options = false;
Append_To_Scatterplot_Title = '';

% Run dataplot and scatplot scripts

[saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Working_Dir, Manuals_Dir);
scatplot(saved_data, drange, ...
         'Manuals_Dir', Manuals_Dir, ...
         'Scatterplot_Inputs_File', Scatterplot_Inputs_File, ...
         'Stats_Output', Stats_Output, ...
         'Output_File', Output_File, ...
         'Statistics_Tex_Output', Statistics_Tex_Output, ...
         'NRC_Options', NRC_Options, ...
         'Append_To_Scatterplot_Title', Append_To_Scatterplot_Title)

Special cases

disp('turb_model...');                    turb_model
disp('jet_decay...');                     jet_decay
disp('wall_model...');                    wall_model
disp('pyrolysis...');                     pyrolysis
disp('birch_tga...');                     birch_tga
disp('water_ice_water...');               water_ice_water
disp('pcm_slab...');                      pcm_slab
disp('pulsating...');                     pulsating
disp('compression_wave...');              compression_wave
disp('plate_view_factor...');             plate_view_factor
disp('terminal_velocity_convergence...'); terminal_velocity_convergence
disp('fluid_part...');                    fluid_part
disp('extinction...');                    extinction
disp('fan_curve...');                     fan_curve
disp('mesh_transformation...');           mesh_transformation
disp('synthetic_eddy_method...');         synthetic_eddy_method
disp('shunn_mms_error...');               shunn_mms_error
disp('openmp_timing_benchmarks...');      openmp_timing_benchmarks
disp('rms_cov_corr...');                  rms_cov_corr
disp('hot_layer_collapse...');            hot_layer_collapse
disp('radiating_polygon...');             radiating_polygon
disp('saad_mms_temporal_error...');       saad_mms_temporal_error
disp('shunn_mms_temporal_error...');      shunn_mms_temporal_error
disp('scaling_tests...');                 scaling_tests
disp('hvac_mass_transport...');           hvac_mass_transport
disp('ht3d_cond...');                     ht3d_cond

display('verification scripts completed successfully!')
