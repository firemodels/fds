% Hostikka
% April 23, 2010
% plate_view_factor.m

close all
clear all

addpath('../../Verification/Radiation')

% set the plot style parameters

plot_style

% Collect data

NRA = [30 60 100];

Exact_Flux_2D = 105.34;
Exact_Flux_cart = 81.8;
Exact_Flux_cyl = 74.1;

% 2D
M = csvread('plate_view_factor_2D_30_devc.csv',2,0);
Flux_2D(1) = max(M(:,2));
M = csvread('plate_view_factor_2D_60_devc.csv',2,0);
Flux_2D(2) = max(M(:,2));
M = csvread('plate_view_factor_2D_100_devc.csv',2,0);
Flux_2D(3) = max(M(:,2));

% Cart
M = csvread('plate_view_factor_cart_30_devc.csv',2,0);
Flux_cart(1) = max(M(:,2));
M = csvread('plate_view_factor_cart_60_devc.csv',2,0);
Flux_cart(2) = max(M(:,2));
M = csvread('plate_view_factor_cart_100_devc.csv',2,0);
Flux_cart(3) = max(M(:,2));

% Cylindrical
M = csvread('plate_view_factor_cyl_30_devc.csv',2,0);
Flux_cyl(1) = max(M(:,2));
M = csvread('plate_view_factor_cyl_60_devc.csv',2,0);
Flux_cyl(2) = max(M(:,2));
M = csvread('plate_view_factor_cyl_100_devc.csv',2,0);
Flux_cyl(3) = max(M(:,2));

% IBM
M = csvread('plate_view_factor_ibm_30_devc.csv',2,0);
Flux_ibm(1) = max(M(:,2));
M = csvread('plate_view_factor_ibm_60_devc.csv',2,0);
Flux_ibm(2) = max(M(:,2));
M = csvread('plate_view_factor_ibm_100_devc.csv',2,0);
Flux_ibm(3) = max(M(:,2));

% Plot data

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y 1.3*Plot_Width Plot_Height])

h=plot(NRA,Exact_Flux_2D*[1 1 1],'r-',NRA,Flux_2D,'ro');
hold on
h=plot(NRA,Exact_Flux_cart*[1 1 1],'b-',NRA,Flux_cart,'bs');
%h=plot(NRA,Flux_ibm,'bd'); % Hold for FDS 7
h=plot(NRA,Exact_Flux_cyl*[1 1 1],'g-',NRA,Flux_cyl,'gd');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

axis([20 110 30 120])
xlabel('Number of Radiation Angles','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Heat Flux (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
Plot_Title='Radiative Heat Flux (plate\_view\_factor)';
text(30,113,Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
%leg=legend('Exact 2D','FDS 2D','Exact 3D','FDS 3D','FDS 3D IBM','Exact cyl.','FDS cyl.','Location','EastOutside');
leg=legend('Exact 2D','FDS 2D','Exact 3D','FDS 3D','Exact cyl.','FDS cyl.','Location','EastOutside');
set(leg,'FontSize',Key_Font_Size);

% add version string if file is available

Git_Filename = 'plate_view_factor_2D_30_git.txt';
addverstr(gca,Git_Filename,'linear')

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[1.3*Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 1.3*Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/plate_view_factor')

%close
