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

%layer_height

Dataplot_Inputs_File = [pwd,'/FDS_validation_dataplot_inputs.csv'];
Working_Dir = [pwd, '/../../Validation/'];
Manuals_Dir = [pwd, '/../../Manuals/'];

Scatterplot_Inputs_File = [pwd, '/FDS_validation_scatterplot_inputs.csv'];
Stats_Output = 'Validation';
Output_File = [pwd, '/FDS_validation_scatterplot_output.csv'];
Statistics_Tex_Output = [pwd, '/../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_statistics.tex'];
Histogram_Tex_Output = [pwd, '/../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_histograms.tex'];
NRC_Options = false;
Append_To_Scatterplot_Title = '';

% The relevant rows from FDS_validation_dataplot_inputs.csv are the last arguments of the following call:

[saved_data,drange] = dataplot(Dataplot_Inputs_File, Working_Dir, Manuals_Dir,'Hamins Gas Burners');
     
scatplot(saved_data, drange, ...
         'Manuals_Dir', Manuals_Dir, ...
         'Scatterplot_Inputs_File', Scatterplot_Inputs_File, ...
         'Stats_Output', Stats_Output, ...
         'Output_File', Output_File, ...
         'Statistics_Tex_Output', Statistics_Tex_Output, ...
         'Histogram_Tex_Output', Histogram_Tex_Output, ...
         'NRC_Options', NRC_Options, ...
         'Append_To_Scatterplot_Title', Append_To_Scatterplot_Title)

display('Hamins_Plots completed successfully!')
