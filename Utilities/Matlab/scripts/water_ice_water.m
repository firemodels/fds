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
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

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

h=plot(wiw_M(1:500,3),wiw_m(1:500,5),'b-',wiw_M(501:1000,3),wiw_m(501:1000,5),'r-');

% Plot attributes
    
set(gca,'FontName',Font_Name)
xlabel('Temperature ($^\circ$C)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Liquid concentration (kg/m$^3$)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
legend('Cooling (freezing)','Heating (melting)')
axis([-10 10 0 1000])
set(h,'LineStyle','-')

% add SVN if file is available

SVN_Filename = ['water_ice_water_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% Create the PDF files

set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/water_ice_water')
