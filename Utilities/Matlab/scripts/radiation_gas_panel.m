
% radiation_gas_panel.m

close all
clear all

% Define standard plotting parameters
plot_style

% Directories
data_dir  = '../../Verification/Radiation/';
plot_dir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% Read the output data
if(~exist([data_dir,'radiation_gas_panel_10cm_devc.csv'])        || ...
   ~exist([data_dir,'radiation_gas_panel_15cm_devc.csv'])        || ...
   ~exist([data_dir,'radiation_gas_panel_25cm_devc.csv'])        || ...
   ~exist([data_dir,'radiation_gas_panel_38cm_devc.csv'])        || ...
   ~exist([data_dir,'radiation_gas_panel_61cm_devc.csv'])        || ...
   ~exist([data_dir,'radiation_gas_panel_76cm_devc.csv'])        || ...
   ~exist([data_dir,'radiation_gas_panel_10cm_offset_devc.csv']) || ...
   ~exist([data_dir,'radiation_gas_panel_15cm_offset_devc.csv']) || ...
   ~exist([data_dir,'radiation_gas_panel_25cm_offset_devc.csv']) || ...
   ~exist([data_dir,'radiation_gas_panel_38cm_offset_devc.csv']) || ...
   ~exist([data_dir,'radiation_gas_panel_61cm_offset_devc.csv']) || ...
   ~exist([data_dir,'radiation_gas_panel_76cm_offset_devc.csv'])) 
    display(['Error One or more of file radiation_gas_panel_<?>_devc.csv does not exist.  Skipping case.']);
    return;
end

M10 = csvread([data_dir,'radiation_gas_panel_10cm_devc.csv'],2);
M15 = csvread([data_dir,'radiation_gas_panel_15cm_devc.csv'],2);
M25 = csvread([data_dir,'radiation_gas_panel_25cm_devc.csv'],2);
M38 = csvread([data_dir,'radiation_gas_panel_38cm_devc.csv'],2);
M61 = csvread([data_dir,'radiation_gas_panel_61cm_devc.csv'],2);
M76 = csvread([data_dir,'radiation_gas_panel_76cm_devc.csv'],2);

M10Offset = csvread([data_dir,'radiation_gas_panel_10cm_offset_devc.csv'],2);
M15Offset = csvread([data_dir,'radiation_gas_panel_15cm_offset_devc.csv'],2);
M25Offset = csvread([data_dir,'radiation_gas_panel_25cm_offset_devc.csv'],2);
M38Offset = csvread([data_dir,'radiation_gas_panel_38cm_offset_devc.csv'],2);
M61Offset = csvread([data_dir,'radiation_gas_panel_61cm_offset_devc.csv'],2);
M76Offset = csvread([data_dir,'radiation_gas_panel_76cm_offset_devc.csv'],2);

% Collect heat fluxes computed by FDS.
xFDS       = [.10,  .15, .25, .38, .46, .61, .76];
yFDS       = [M10(end,end),       M15(end,end),       M25(end,end),       M38(end,end),       M46(end,end),       M61(end,end),       M76(end,end)];
yFDSOffset = [M10Offset(end,end), M15Offset(end,end), M25Offset(end,end), M38Offset(end,end), M46Offset(end,end), M61Offset(end,end), M76Offset(end,end)];

% On- and off-axis heat fluxes computed using Boltzmann's law and configuration factors.
xCF       = [.10,  .15,  .25,  .38,  .46,  .61,  .76 ];
yCF       = [71.7, 54.4, 30.8, 16.3, 11.8, 7.10, 4.70];
yCFOffset = [38.4, 24.0, 17.4, 11.5, 9.0,  5.97, 4.18];

% Heat fluxes reported by Simms.
xSimms = [.10,  .15,  .25,  .38,  .46,  .61,  .76 ];
ySimms = [70.2, 51.8, 31.5, 16.9, 12.7, 7.61, 5.52];

% Plot heat flux versus distance.
fig=figure();
ax1=get(fig,'CurrentAxes');

plot(xSimms,ySimms,    '-o', 'MarkerSize', 10, 'MarkerFaceColor', 'red'); hold on
plot(xCF,   yCF,       '-o', 'MarkerSize', 10, 'MarkerFaceColor', 'green'); hold on
plot(xFDS,  yFDS,      '-o', 'MarkerSize', 10, 'MarkerFaceColor', 'blue'); hold on

xlabel('Perpendicular distance from gas panel (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name); hold on
ylabel('Heat flux (kW/m^2)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlim([0, 0.8])
ylim([0, 80])

xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8])

lh1=legend({'Simms', 'Configuration factor', 'FDS'}, 'Position', [0.6 0.7 0.3 0.2]);
lh1.Title.String = 'On-axis';
set(lh1,'FontSize',Key_Font_Size)

% Matlab will draw only one legend per set of axes.  
% Draw a second set of axes `under' the first to
% overcome this limitation.
ax1=gca;

ax2=axes('Position',get(ax1,'Position'),'xlim',get(ax1,'xlim'),'ylim',get(ax1,'ylim'),'Visible','off','Color','none');

hold(ax2,'on'); % "do not erase" when plotting

set(ax2,'FontName',Font_Name)
set(ax2,'FontSize',Title_Font_Size)

plot(ax2, xCF,   yCFOffset,  '-s', 'MarkerSize', 10, 'MarkerFaceColor', 'green');
plot(ax2, xFDS,  yFDSOffset, '-s', 'MarkerSize', 10, 'MarkerFaceColor', 'blue');

lh2=legend(ax2, {'Configuration factor', 'FDS'}, 'Position', [0.6 0.5 0.3 0.15]);
lh2.Title.String='Off-axis';
set(lh2,'FontSize',Key_Font_Size)

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
