% Hostikka
% 11-15-2010
% water_ice_water.m
%
% This reads in and plots the data for water_ice_water test case.

close all
clear all

% Set global reaction rate parameters

addpath('../../../Verification/Pyrolysis')

close all
    
plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

wiw_M = csvread('water_ice_water_devc.csv',2);
wiw_m = csvread('water_ice_water_prof_01.csv',2);

h=plot(wiw_M(1:500,3),wiw_m(1:500,5),wiw_M(501:1000,3),wiw_m(501:1000,5),'r-');

% Plot attributes
    
set(gca,'FontName',Font_Name)
xlabel('Temperature ($^\circ$C)','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('Liquid concentration (kg/m$^3$)','Interpreter','LaTeX','FontSize',Label_Font_Size)
legend('Cooling','Heating')
axis([-10 10 0 1000])
set(h,'LineStyle','-')

% Create the PDF files

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/water_ice_water')
