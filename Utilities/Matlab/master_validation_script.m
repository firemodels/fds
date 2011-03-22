% McGrattan
% March 21, 2011
% master_validation_script.m

% This script creates the plots that are included in the FDS Validation 
% Guide. It consists of calls to other scripts contained within the
% subfolder called "scripts". 

% The most important script is called dataplot. It reads the file called
% validation_data_config_matlab.csv and generates 1000+ plots. If you
% want to process only some of these plots, comment out the other 
% scripts and change the data plot line as follows:

% [saved_data,drange] = dataplot('validation',[a:b]);

% where a and b are the lines in the .csv file you want to process.

close all
clear all
addpath 'scripts'

% Scripts that run prior to dataplot

flame_height

% dataplot creates most of the plots for the Validation Guide. It must
% be run before scatplot, which makes the scatter plots.

[saved_data,drange] = dataplot('validation');
scatplot(saved_data,drange)

% Miscellaneous other scripts for special cases

beyler_hood
check_hrr
sandia_helium_plume
sandia_methane_fire
BRE_spray
% Cup_burner
vettori_flat
vettori_sloped
 
display('validation scripts completed successfully!')
