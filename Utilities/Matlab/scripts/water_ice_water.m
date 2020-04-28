% Hostikka
% 11-15-2010
% water_ice_water.m
%
% This reads in and plots the data for water_ice_water test case.

close all
clear all

% Set global reaction rate parameters

addpath('../../Verification/Pyrolysis')

close all

plot_style

skip_case = 0;
if ~exist('water_ice_water_devc.csv')
    display('Error: File water_ice_water_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('water_ice_water_prof_01.csv')
    display('Error: File water_ice_water_prof_01.csv does not exist. Skipping case.')
    skip_case = 1;
end
if skip_case
    return
end

wiw_M = csvread('water_ice_water_devc.csv',2);
wiw_m = csvread('water_ice_water_prof_01.csv',2);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

h=plot(wiw_M(1:500,3),wiw_m(1:500,5),'b-',wiw_M(501:1000,3),wiw_m(501:1000,5),'r-');

% Plot attributes

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Liquid concentration (kg/m^3)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
axis([-10 10 0 1000])

lh=legend('Cooling (freezing)','Heating (melting)');
set(lh,'FontSize',Key_Font_Size,'Interpreter',Font_Interpreter)

% add git version if file is available

git_file = ['water_ice_water_git.txt'];
addverstr(gca,git_file,'linear')

% Create the PDF files

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/water_ice_water')
