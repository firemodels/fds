% McGrattan
% March 21, 2011
% FDS_validation_script.m
%
% This script creates the plots that are included in the FDS Validation 
% Guide. It consists of calls to other scripts contained within the
% subfolder called "scripts". 
%
% The most important script is called dataplot. It reads the file called
% FDS_validation_dataplot_inputs.csv and generates 1000+ plots. If you
% want to process only some of these plots, comment out the other 
% scripts and change the data plot line as follows:
%
% [saved_data,drange] = dataplot(Dataplot_Inputs_File,Validation_Dir,Manuals_Dir,[a:b]);
%
% where a and b are the lines in the .csv file you want to process.
% Alternatively, you can specify the "Dataname" you want:
%
% [saved_data,drange] = dataplot(Dataplot_Inputs_File,Validation_Dir,Manuals_Dir,'WTC');
%
% In this case, all the WTC results will be plotted.
%
% Dataplot creates most of the plots for the Validation Guide.
% It must be run before scatplot, which makes the scatter plots.

close all
clear all

addpath 'scripts'

% Scripts that run prior to dataplot

flame_height
cat_mccaffrey
NIST_RSE
sippola_aerosol_deposition
layer_height
Hamins_CH4_Average
combine_csiro
fm_datacenter_scatter

% Dataplot and scatplot options

Dataplot_Inputs_File = [pwd,'/FDS_validation_dataplot_inputs.csv'];
Working_Dir = [pwd, '/../../Validation/'];
Manuals_Dir = [pwd, '/../../FDS/Manuals/'];
Scatterplot_Inputs_File = [pwd, '/FDS_validation_scatterplot_inputs.csv'];

% Statistics output options

Stats_Output = 'Validation';
Output_File = [pwd, '/FDS_validation_scatterplot_output.csv'];
Statistics_Tex_Output = [pwd, '/../../FDS/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_statistics.tex'];
Histogram_Tex_Output = [pwd, '/../../FDS/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/ScatterPlots/validation_histograms.tex'];

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
         'Histogram_Tex_Output', Histogram_Tex_Output, ...
         'NRC_Options', NRC_Options, ...
         'Append_To_Scatterplot_Title', Append_To_Scatterplot_Title)
     
% Miscellaneous other scripts for special cases

backward_facing_step
beyler_hood
sandia_helium_plume
sandia_methane_fire
spray_attenuation
Cup_burner
flame_height2
purdue_flames
christifire
pressure_coefficient
VTT_Sprays
fm_datacenter_veltest
umd_line_burner

display('validation scripts completed successfully!')
