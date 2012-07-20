%-------------------------------------------------------------------------%
% FILE: sawtooth.m
% AUTHOR: Max Gould
% LAST DATE EDITED: 19 July 2012
%-------------------------------------------------------------------------%

%-------------------------------------------------------------------------%
% Generate velocity plot for SAWTOOTH=.FALSE. example in the FDS User Guide
%-------------------------------------------------------------------------%
close all
clear all

% Open directory with input files
input_dir = '../../Verification/Flowfields/';

%-------------------------------------------------------------------------%
% Define Parameters and Variables
%-------------------------------------------------------------------------%

% Simulation time
TIME = 5.0;
ROW = 100*TIME+2;

% Plot range
Index = 0:31;

%-------------------------------------------------------------------------%
% Import Data
%-------------------------------------------------------------------------%

SawtoothFIn = csvread([input_dir,'sawtooth_false_devc.csv'],ROW,1,...
                      [ROW,1,ROW,32]);
SawtoothTIn  = csvread([input_dir,'sawtooth_true_devc.csv'],ROW,1,...
                      [ROW,1,ROW,32]);

%-------------------------------------------------------------------------%
% Plot Data
%-------------------------------------------------------------------------%

% Open directory for plot output
plot_dir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';

% Define standard plotting parameters
plot_style

% Produce plots
SawtoothFPlot = plot(Index,SawtoothFIn);
hold on
SawtoothTPlot = plot(Index,SawtoothTIn);
hold off

% Define plot colors
set(SawtoothFPlot,'Color',[1.0 0.8 0.0]);
set(SawtoothTPlot,'Color',[0.0 0.8 1.0]);

% Specify axes for above plots
axis([0.0 31.0 0.0 2.4]);

% Set standard plotting parameters
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);

% Label axes
xlabel('Z-Coordinate (mm)','FontSize',Title_Font_Size,...
       'Interpreter',Font_Interpreter);
ylabel('U-Velocity (m/s)','FontSize',Title_Font_Size,...
       'Interpreter',Font_Interpreter);

% Create legend
leg1 = legend('Sawtooth = False','Sawtooth = True','Location','NorthEast');
set(leg1,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter);

% Save plot to file
print(gcf,'-dpdf',[plot_dir,'sawtooth_ugraph']);



