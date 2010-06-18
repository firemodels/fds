% Hostikka
% April 23, 2010
% plate_view_factor.m

close all
clear all

addpath('../../../Verification/Radiation')

plot_style

% Collect data

NRA = [30 60 100];

Exact_Flux_2D = 105.34;
Exact_Flux_cart = 81.8;
Exact_Flux_cyl = 74.1; 

%2D 
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
%cylindrical
M = csvread('plate_view_factor_cyl_30_devc.csv',2,0);
Flux_cyl(1) = max(M(:,2));
M = csvread('plate_view_factor_cyl_60_devc.csv',2,0);
Flux_cyl(2) = max(M(:,2));
M = csvread('plate_view_factor_cyl_100_devc.csv',2,0);
Flux_cyl(3) = max(M(:,2));

% Plot data

%2D

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)

h=plot(NRA,Exact_Flux_2D*[1 1 1],'r-',NRA,Flux_2D,'ro--');
hold on
h=plot(NRA,Exact_Flux_cart*[1 1 1],'b-',NRA,Flux_cart,'bo--');
h=plot(NRA,Exact_Flux_cyl*[1 1 1],'g-',NRA,Flux_cyl,'go--');

axis([20 110 40 120])
set(gca,'FontName',Font_Name)
xlabel('Number radiation angles','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('Radiative heat flux (kW/m$^2$)','Interpreter','LaTeX','FontSize',Label_Font_Size)
title('Radiative heat flux (plate\_view\_factor)','Interpreter','LaTeX','FontSize',Label_Font_Size)
legend('Exact 2D','FDS 2D','Exact cart.','FDS cart.','Exact cyl.','FDS cyl.',...
   'Location','SouthEast')

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_5_Verification_Guide/SCRIPT_FIGURES/plate_view_factor')

%close