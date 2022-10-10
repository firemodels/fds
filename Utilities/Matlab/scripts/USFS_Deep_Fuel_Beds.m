
% 6/24/2022 Noelle Crump

% This code is set up to automate the production of a set of
% simulation Validation plots.  A list of parameters for each experiment is read
% in from a file. The units are not standard SI.

% Here are the columns:
% 1 run# -- a sequential number to track the simulation.
% 2 burn# -- the number of the actual deep fuel beds burn, if applicable
% 3 spacing -- between center points of fuel arrays in cm.
% 4 angle -- angle in degrees of inclined platform.
% 5 depth -- in inches.
% 6 moisture content -- in percent, dry weight basis, so mc = 5 = (wet-dry)/(dry).
% 9 temp -- ambient temp in degrees fahrenheit
% 10 rh -- relative humidity, in percent
% 11 burned -- 1 or zero; 1 means "burned", zero means "did not burn".
% 12 pct_burn -- percent of fuel burned.

close all
clear all

% Runs from Matlab folder (Validation plots implementation)
outdir = '../../../out/USFS_Deep_Fuel_Beds/';
expdir = '../../../exp/USFS_Deep_Fuel_Beds/';
plot_dir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/USFS_Deep_Fuel_Beds/';

EXP = importdata([expdir,'exp_params.csv'],',',1);

BURN_NO      = EXP.data(:,find(strcmp(EXP.colheaders,'BURN_NO')));
SPACING      = EXP.data(:,find(strcmp(EXP.colheaders,'SPACING'))).*2.5;
DEPTH        = EXP.data(:,find(strcmp(EXP.colheaders,'DEPTH'))).*2.5;
FMC          = EXP.data(:,find(strcmp(EXP.colheaders,'MOISTURE_CONTENT')));
SLOPE        = EXP.data(:,find(strcmp(EXP.colheaders,'SLOPE')));
BURN_PCT_EXP = EXP.data(:,find(strcmp(EXP.colheaders,'BURN_PERCENT')));

% read in data
for i = 1:length(BURN_NO)

    CHID = ['burn',num2str(BURN_NO(i)),'_',num2str(DEPTH(i)),'D_',num2str(SLOPE(i)),'S_',num2str(SPACING(i)),'L'];
    HRR  = importdata([outdir,[CHID,'_cat_hrr.csv']],',', 2);
    DEVC = importdata([outdir,[CHID,'_cat_devc.csv']],',', 2);

    % Determine FDS (out) Burn Type:
    hrr_col = HRR.data(:,(strcmp(HRR.colheaders,'HRR'))); % (kw)
    time_col = HRR.data(:,(strcmp(HRR.colheaders,'Time'))); % (s)
    fuel_load = DEVC.data(2,(strcmp(DEVC.colheaders,'"Fuel Load"'))); % (kg)
    dry_fuel_mass = (1/(1+FMC(i)*0.01))*fuel_load; % (kg)
    total_fuel_energy = dry_fuel_mass*17425.0; % (kJ) 17425 = heat of combustion

    ignition_energy = trapz([0,10,20],[500*0.3*1.6,500*0.3*1.6,0]); % (kJ) hrr ramp
    total_fds_energy = trapz (time_col,hrr_col); % (kJ)

    BURN_PCT_FDS(i) = (total_fds_energy-ignition_energy)/total_fuel_energy * 100.;
end

figure
plot_style

set(gca,'Units',Plot_Units)
set(gca,'Position',[Scat_Plot_X Scat_Plot_Y Scat_Plot_Width Scat_Plot_Height])

h = scatter(BURN_PCT_EXP,BURN_PCT_FDS,10,'filled');
box on

hold on
plot([0,100],[0,100],'k-')

for i=1:length(BURN_NO)
    % add run number to each marker
    text(BURN_PCT_EXP(i),BURN_PCT_FDS(i),['     ',num2str(BURN_NO(i))],'FontName',Font_Name,'FontSize',2);
end

xlim([0,100])
ylim([0,100])

xlabel('Fuel Energy Burned EXP (%)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)
ylabel('Fuel Energy Burned FDS (%)','Interpreter',Font_Interpreter,'FontSize',Scat_Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Scat_Label_Font_Size)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'Position',[0 0 Scat_Paper_Width Scat_Paper_Height]);
print(gcf,'-dpdf',[plot_dir,'USFS_Deep_Fuel_Beds_Scatterplot']);

return