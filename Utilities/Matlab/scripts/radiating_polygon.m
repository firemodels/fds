% McGrattan
% 2-16-15
% radiating_polygon.m

close all
clear all

datadir='../../Verification/Radiation/';

skip_case=false;

if ~exist([datadir,'radiating_polygon_square_20_line.csv']); skip_case=true; end
if ~exist([datadir,'radiating_polygon_square_40_line.csv']); skip_case=true; end
if ~exist([datadir,'radiating_polygon_square_80_line.csv']); skip_case=true; end

if skip_case % skip_case_if
    display(['Error: Files for radiating_polygon_square do not exist. Skipping case.'])
else

% gather FDS results

M = importdata([datadir,'radiating_polygon_square_20_line.csv'],',',2); z_20 = 1-M.data(:,1); flux_20 = M.data(:,2);
M = importdata([datadir,'radiating_polygon_square_40_line.csv'],',',2); z_40 = 1-M.data(:,1); flux_40 = M.data(:,2);
M = importdata([datadir,'radiating_polygon_square_80_line.csv'],',',2); z_80 = 1-M.data(:,1); flux_80 = M.data(:,2);

% analytical solution (Siegel and Howell, 2nd ed., appendix)

n=4;
pi=4.*atan(1);

for j=1:1000
   z(j)=1-0.001*(j-1);
   R=0.5*sqrt(2)/z(j);
   H=0.5/z(j);
   flux(j)=5.67e-11*1273.15^4*n*H/(pi*sqrt(1+H^2))*atan(sqrt((R^2-H^2)/(1+H^2)));
end

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
Marker_Size=1;

H(1)=plot(z,flux,'k-'); hold on
H(2)=plot(z_20,flux_20,'-mo','MarkerSize',Marker_Size);
H(3)=plot(z_40,flux_40,'-go','MarkerSize',Marker_Size);
H(4)=plot(z_80,flux_80,'-ro','MarkerSize',Marker_Size);

axis([0 1 30 150])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

text(0.25,140,'Radiative Flux from a Hot Square Plate','FontSize',Label_Font_Size,'FontName',Font_Name)

xlabel('Distance from Plate (m)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name);
ylabel('Heat Flux (kW/m^2)','FontSize',Label_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
legend_handle = legend(H,'analytical','5 cm, 100 angles','2.5 cm, 400 angles','1.25 cm, 1600 angles','Location','SouthWest');
set(legend_handle,'Interpreter',Font_Interpreter);
set(legend_handle,'Fontsize',Label_Font_Size);
set(legend_handle,'Box','on');

% add Git revision if file is available

Git_Filename = [datadir,'radiating_polygon_square_20_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/radiating_polygon_square'])

end % skip_case_if
