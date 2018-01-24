%
% This test is now being driven from dataplot.  Retain
% the file as an example
%

close all
clear all

% Define standard plotting parameters
plot_style

% Directories
data_dir  = '../../Verification/Vegetation/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Read the output data
if ~exist([data_dir,'vege_mass_conservation_devc.csv']) 
    display(['Error File vege_mass_conservation_devc.csv does not exist.  Skipping case.']);
    return;
end

M = csvread([data_dir,'vege_mass_conservation_devc.csv'],2);

% Total mass should sum to 1.0 at every time
mass_delta = 1.0 - sum( M(:,[2:4]), 2 );

% Plot mass_delta versus time, M(:,1)
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(M(:,1),mass_delta,'LineWidth',Line_Width); hold on

xlabel('Time (s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
ylabel('\Delta mass from expected (kg)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name)
axis([0 max(M(:,1)) -0.5 1])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% Add Git revision if file is available

Git_Filename = [data_dir,'vege_mass_conservation_git.txt'];
addverstr(gca,Git_Filename,'linear')

hold off

% Print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'vege_mass_conservation'])
