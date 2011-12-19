% McGrattan
% March 21, 2011
% master_validation_script.m
%
% This script creates the plots that are included in the FDS Validation 
% Guide. It consists of calls to other scripts contained within the
% subfolder called "scripts". 
%
% The most important script is called dataplot. It reads the file called
% validation_data_config_matlab.csv and generates 1000+ plots. If you
% want to process only some of these plots, comment out the other 
% scripts and change the data plot line as follows:
%
% [saved_data,drange] = dataplot(cfil,vdir,plotdir,[a:b]);
%
% where a and b are the lines in the .csv file you want to process.
% Alternatively, you can specify the "Dataname" you want:
%
% [saved_data,drange] = dataplot(cfil,vdir,plotdir,'WTC');
%
% In this case, all the WTC results will be plotted. 

close all
clear all

addpath 'scripts'

% Name of dataplot input file, directory where data is located, and 
% directory to put the plots

cfil = [pwd,'/validation_data_config_matlab.csv'];
vdir = [pwd,'/../../Validation/'];
plotdir = [pwd,'/../../Manuals/'];

% Scripts that run prior to dataplot

flame_height
cat_mccaffrey
NIST_RSE

% dataplot creates most of the plots for the Validation Guide. It must be run before scatplot, which makes the scatter plots.

[saved_data,drange] = dataplot(cfil,vdir,plotdir);
scatplot(saved_data,drange)

% Miscellaneous other scripts for special cases

harrisonplumes
beyler_hood
check_hrr
sandia_helium_plume
sandia_methane_fire
spray_attenuation
Cup_burner
vettori_flat
vettori_sloped
flame_height2
 
display('validation scripts completed successfully!')
