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

restoredefaultpath
addpath 'scripts'

% Scripts to run prior to dataplot

disp('radiation_plane_layer...');          radiation_plane_layer
disp('ns2d...');                           ns2d
disp('vort2d...');                         vort2d
disp('flame_species...');                  flame_species
disp('convective_cooling_convergence...'); convective_cooling_convergence
disp('random_walk_soln...');               random_walk_soln
disp('vegetation_absorb...');              vegetation_absorb

% Dataplot and scatplot options

Dataplot_Inputs_File = 'FDS_verification_dataplot_inputs.csv';
Working_Dir = '../../Verification/';
Manuals_Dir = '../../Manuals/';
Scatterplot_Inputs_File = 'FDS_verification_scatterplot_inputs.csv';

% Statistics output options

Stats_Output = 'Verification';
Scatterplot_Dir = [pwd, '/../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/Scatterplots/'];

% Run dataplot and scatplot scripts

[saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Working_Dir, Manuals_Dir);
scatplot(saved_data, drange, ...
         'Manuals_Dir', Manuals_Dir, ...
         'Scatterplot_Inputs_File', Scatterplot_Inputs_File, ...
         'Stats_Output', Stats_Output, ...
         'Scatterplot_Dir', Scatterplot_Dir)

% Special cases

disp('compression_wave...');              compression_wave
disp('extinction...');                    extinction
disp('rotcube_cc_mms_error...');          rotcube_cc_mms_error
disp('openmp_timing_benchmarks...');      openmp_timing_benchmarks
disp('rms_cov_corr...');                  rms_cov_corr
disp('hot_layer_collapse...');            hot_layer_collapse
disp('radiating_polygon...');             radiating_polygon
disp('favre_test...');                    favre_test
disp('scaling_tests...');                 scaling_tests
disp('particle_size_distribution...');    particle_size_distribution
disp('mass_balance_reac...');             mass_balance_reac
disp('mass_balance_gas_volume...');       mass_balance_gas_volume
disp('ht3d_sphere...');                   ht3d_sphere
disp('geom_positive_errors...');          geom_positive_errors
disp('geom_channel_test...');             geom_channel_test
disp('atmospheric_boundary_layer...');    atmospheric_boundary_layer
disp('htc_forced...');                    htc_forced
disp('natconh...');                       natconh
disp('natconv...');                       natconv
disp('freecon_sphere...');                freecon_sphere
disp('nat_conv_hot_plate...');            nat_conv_hot_plate
disp('impinging_jet...');                 impinging_jet

display('verification scripts completed successfully!')
