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

skip_case = 0;

if ~exist('plate_view_factor_2D_30_devc.csv')
    display('Error: File plate_view_factor_2D_30_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('plate_view_factor_2D_60_devc.csv')
    display(['Error: File plate_view_factor_2D_60_devc.csv does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist('plate_view_factor_2D_100_devc.csv')
    display(['Error: File plate_view_factor_2D_100_devc.csv does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist('plate_view_factor_cart_30_devc.csv')
    display('Error: File plate_view_factor_cart_30_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('plate_view_factor_cart_60_devc.csv')
    display(['Error: File plate_view_factor_cart_60_devc.csv does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist('plate_view_factor_cart_100_devc.csv')
    display(['Error: File plate_view_factor_cart_100_devc.csv does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist('plate_view_factor_cyl_30_devc.csv')
    display('Error: File plate_view_factor_cyl_30_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('plate_view_factor_cyl_60_devc.csv')
    display(['Error: File plate_view_factor_cyl_60_devc.csv does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist('plate_view_factor_cyl_100_devc.csv')
    display(['Error: File plate_view_factor_cyl_100_devc.csv does not exist. Skipping case.'])
    skip_case = 1;
end


if skip_case
    return
end

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

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

h=plot(NRA,Exact_Flux_2D*[1 1 1],'r-',NRA,Flux_2D,'ro');
hold on
h=plot(NRA,Exact_Flux_cart*[1 1 1],'b-',NRA,Flux_cart,'bs');
h=plot(NRA,Exact_Flux_cyl*[1 1 1],'g-',NRA,Flux_cyl,'gd');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

axis([20 110 30 120])
xlabel('Number radiation angles','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Radiative heat flux (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
Plot_Title='Radiative heat flux (plate\_view\_factor)';
text(30,113,Plot_Title,'FontSize',Title_Font_Size,'FontName',Font_Name,'Interpreter',Font_Interpreter)
leg=legend('Exact 2D','FDS 2D','Exact cart.','FDS cart.','Exact cyl.','FDS cyl.',...
   'Location','SouthEast');
set(leg,'FontSize',Key_Font_Size);

% add version string if file is available

Git_Filename = 'plate_view_factor_2D_30_git.txt';
addverstr(gca,Git_Filename,'linear')

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/plate_view_factor')

%close