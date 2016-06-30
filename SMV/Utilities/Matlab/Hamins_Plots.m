% McGrattan
% June 14, 2016
% Hamins_Plots.m
%
% This script creates the plots of Hamins' burner experiments.
% 
% The data is in the FDS-SMV repository, Validation/Hamins_CxHy/Experimental_Data
%
% The plots are put in Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Hamins_CxHy
%
% The plot details are in Utilities/Matlab/FDS_validation_dataplot_inputs.csv starting at row 2664
%

close all
clear all

addpath 'scripts'

Dataplot_Inputs_File = [pwd,'/FDS_validation_dataplot_inputs.csv'];
Working_Dir = [pwd, '/../../Validation/'];
Manuals_Dir = [pwd, '/../../Manuals/'];

% The relevant rows from FDS_validation_dataplot_inputs.csv are the last arguments of the following call:

[saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Manuals_Dir,[2664:2828]);
     
display('Hamins_Plots completed successfully!')
