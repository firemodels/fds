
% radiation_gas_panel.m

close all
clear all

% Define standard plotting parameters
plot_style

% Directories
data_dir  = '../../Verification/Radiation/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Read the output data
if ~exist([data_dir,'radiation_gas_panel_devc.csv']) 
    display(['Error File radiation_gas_panel_devc.csv does not exist.  Skipping case.']);
    return;
end

M = csvread([data_dir,'radiation_gas_panel_devc.csv'],2);

% Collect heat fluxes computed by FDS.
xFDS = [.10,  .15, .25, .38, .61, .76];
yFDS = [4,5,6];

% Heat fluxes computed using Boltzmann's law and configuration
% factors.
xCF = [.10,  .15,    .25,  .38,  .46,  .76];
yCF = [71.7, 564.4,  30.8, 16.3, 11.8, 4.70];

% Heat fluxes reported by Simms.
xSimms = [.10,  .15,  .25,  .38,  .61, ..76];
ySimms = [69.4, 53.3, 32.2, 17.0, 7.45, 4.92];

% Plot heat flux versus distance.
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(xFDS,  yFDS,  '-o', 'MarkerSize', 15, 'MarkerFaceColor', 'blue' ); hold on
plot(xCF,   yCF,   '-s', 'MarkerSize', 15, 'MarkerFaceColor', 'green'); hold on
plot(xSimms,ySimms,'-d', 'MarkerSize', 15, 'MarkerFaceColor', 'red' .);

xlabel('Distance from center of gas panel (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
ylabel('Heat flux (kW/m^2)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
axis([0 max(M(:,1)) -0.5 1])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% Add Git revision if file is available

Git_Filename = [data_dir,'radiation_gas_panel_git.txt'];
addverstr(gca,Git_Filename,'linear')

hold off

% Print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'radiation_gas_panel'])
