
% natural_convective_cooling_Frankman.m
% See "Fine fuel heating by radiant flux", Frankman et al., 
% Combustion Science and Technology, 2010, 182:215-230.

close all
clear all

% Define standard plotting parameters
plot_style

% Directories
data_dir  = '../../../Verification/Heat_Transfer/';
plot_dir = '../../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frankman data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The rise in steady state fuel surface temperatures in deg C 
% estimated from Frankman Figure 4.  While the absolute fuel
% temperatures are reported in Table 1, the ambient temperature
% for each experiment is not, hence the need to eyeball Figure 4.
ySE_Frankman_DeltaT = [62., 27., 25., 9. ];
yPP_Frankman_DeltaT = [87., 47., 27., 18.];
yLE_Frankman_DeltaT = [98., 47., 33., 22.];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDS data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read the steady state surface temperature calculated by FDS
% for each of the fine fuels positioned at 15, 25, 35 and 45 cm.
%   SE = Small Excelsior
%   PP = Ponderosa Pine
%   LE = Large Excelsior
if(~exist([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_15cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_25cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_35cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_45cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_15cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_25cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_35cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_45cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_15cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_25cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_35cm_devc.csv']) || ...
   ~exist([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_45cm_devc.csv'])) 
    display(['Error One or more of file Frankman_XptYmmrod_<?>cm_devc.csv does not exist. Skipping case.']);
    return;
end

SE15 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_15cm_devc.csv'],2);
SE25 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_25cm_devc.csv'],2);
SE35 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_35cm_devc.csv'],2);
SE45 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt44mmrod_45cm_devc.csv'],2);

PP15 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_15cm_devc.csv'],2);
PP25 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_25cm_devc.csv'],2);
PP35 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_35cm_devc.csv'],2);
PP45 = csvread([data_dir,'natural_convective_cooling_Frankman_0pt70mmrod_45cm_devc.csv'],2);

LE15 = csvread([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_15cm_devc.csv'],2);
LE25 = csvread([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_25cm_devc.csv'],2);
LE35 = csvread([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_35cm_devc.csv'],2);
LE45 = csvread([data_dir,'natural_convective_cooling_Frankman_1pt29mmrod_45cm_devc.csv'],2);

% The steady state fuel surface temperatures in deg C computed by
% FDS for the fuels positioned at .15, .25, .35 and .45 m from the 
% radiant panel.
x =    [.15, .25, .35, .45];
ySE_FDS  = [SE15(end,end), SE25(end,end), SE35(end,end), SE45(end,end)];
yPP_FDS  = [PP15(end,end), PP25(end,end), PP35(end,end), PP45(end,end)];
yLE_FDS  = [LE15(end,end), LE25(end,end), LE35(end,end), LE45(end,end)];

% Compute the fuel surface temperature *rise* from ambient = 25C
% for the FDS calculations.
ySE_FDS_DeltaT = ySE_FDS - 25.;
yPP_FDS_DeltaT = yPP_FDS - 25.;
yLE_FDS_DeltaT = yLE_FDS - 25.;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot heat flux versus distance.
fig=figure();
ax1=get(fig,'CurrentAxes');

% Generate dummy info to produce the desired plot legend below.
% This a very awkward hack.  The MarkerSize here (7) is
% deliberately set smaller than the MarkerSize displayed in the
% plot (10) in order to hide the symbols generated. The plot
% symbols displayed in the legend are consequently a bit small. 
h = zeros(3,1);
h(1) = plot(x+0.01, ySE_Frankman_DeltaT,'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'white'); hold on;
h(2) = plot(x+0.01, ySE_Frankman_DeltaT,'k^', 'MarkerSize', 7, 'MarkerFaceColor', 'white'); hold on;
h(3) = plot(x+0.01, ySE_Frankman_DeltaT,'ks', 'MarkerSize', 7, 'MarkerFaceColor', 'white'); hold on;

plot(x+0.01, ySE_Frankman_DeltaT, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'green'); hold on;
plot(x+0.01, yPP_Frankman_DeltaT, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'green'); hold on;
plot(x+0.01, yLE_Frankman_DeltaT, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'green'); hold on;

plot(x, ySE_FDS_DeltaT, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'cyan'); hold on;
plot(x, yPP_FDS_DeltaT, 'k^', 'MarkerSize', 10, 'MarkerFaceColor', 'cyan'); hold on;
plot(x, yLE_FDS_DeltaT, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'cyan'); hold on;

% Grid
grid on

% Labels
xlabel('Distance from radiant panel (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name); hold on
ylabel('Fuel temperature rise (deg C)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontName',Font_Name); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

% Limits
xlim([0, 0.5])
ylim([0, 120])

% Ticks
xticks([0.05, 0.15, 0.25, 0.35, 0.45])
yticks([20, 40, 60, 80, 100, 120])

set(gca,'XMinorTick','off','YMinorTick','on')

% Legend
lh1=legend({'Small Excelsior', 'Ponderosa Pine', 'Large Excelsior'}, 'Location', 'northeast');
set(lh1,'FontSize',Key_Font_Size)

text(0.34, 91, 'GREEN =',   'Color', 'green', 'FontWeight', 'bold');
text(0.41, 91, 'Frankman',  'Color', 'black');
text(0.34, 85, 'CYAN   = ', 'Color', 'cyan',  'FontWeight', 'bold');
text(0.41, 85, 'FDS',       'Color', 'black');

% Add Git revision if file is available

Git_Filename = [data_dir,'Frankman_natural_convection_git.txt'];
addverstr(gca,Git_Filename,'linear')

hold off

% Print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'natural_convective_cooling_Frankman'])
