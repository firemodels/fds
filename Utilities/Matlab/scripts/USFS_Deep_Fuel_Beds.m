
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
params_loc = '../../../exp/USFS_Deep_Fuel_Beds/exp_params.csv';
outdir = '../../../out/USFS_Deep_Fuel_Beds/';
expdir = '../../../exp/USFS_Deep_Fuel_Beds/';
plot_dir = [pwd, '/../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/USFS_Deep_Fuel_Beds/'];
     % Runs from scripts folder (testing implementation)
% params_loc = '../../../../exp/USFS_Deep_Fuel_Beds/exp_params.csv';
% outdir = '../../../../out/USFS_Deep_Fuel_Beds/';
% expdir = '../../../../exp/USFS_Deep_Fuel_Beds/';
% plot_dir = [pwd, '/../../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/USFS_Deep_Fuel_Beds/'];

experiment_params = importdata(params_loc, ",",1);
experiment_params = experiment_params.data;
s_experiment_params = size(experiment_params);

no_match_data = [];
match_data = [];
all_percent = [];
    % read in data
for runno = 1:s_experiment_params(1)
    thisline = experiment_params(runno,:);

    %runno = thisline(1); % sequential number for batch file run
    burn_no = thisline(2); % 
    spacing = thisline(3).*0.0254; % in-> m
    angle_r = deg2rad(thisline(4));
    depth = thisline(5).*0.0254; % in -> m
    mstcnt = thisline(6).*.01; % 
    burned = thisline(11); % Did the fire spread? 1 = yes, 0 = no
    pct_burn = thisline(12);
    s_spng = num2str(thisline(3).*2.5);
    s_dpth = num2str(thisline(5).*2.5);
    hrr_filename = ['burn',num2str(burn_no),'_',s_dpth,'D_',num2str(thisline(4)),'S_',s_spng,'L_cat_hrr.csv'];
    devc_filename = ['burn',num2str(burn_no),'_',s_dpth,'D_',num2str(thisline(4)),'S_',s_spng,'L_cat_devc.csv'];
    hrr_data = importdata([outdir,hrr_filename], ",", 2);
    devc_data = importdata([outdir,devc_filename], ",", 2);

    % Do some calculations
    eff_spacing = (spacing * cos(angle_r))*100;      % effective spacing (cm)
    eff_depth = (spacing * sin(angle_r)+depth)*100 ; % effective height of bed (cm)


    % Determine FDS (out) Burn Type:
    hrr_col = hrr_data.data(:,(strcmp(hrr_data.colheaders,'HRR'))); % (kw)
    time_col = hrr_data.data(:,(strcmp(hrr_data.colheaders,'Time'))); % (s)
    fuel_load = devc_data.data(2,(strcmp(devc_data.colheaders,'"Fuel Load"'))); % (kg)
    dry_fuel_mass = (1/(1+mstcnt))*fuel_load; % (kg)
    total_fuel_energy = dry_fuel_mass*17425.0; % (kJ) 17425 = heat of combustion

    ignition_hrr = trapz([0,10,20],[0,(500*.072),0]); % (kJ) hrr ramp
    total_fds_hrr = trapz (time_col,hrr_col); % (kJ)

    percent_consumption = (total_fds_hrr-ignition_hrr)/(total_fuel_energy);

    % use hrr to determine the burn quality
    if percent_consumption > 0.60
        burn_type_o = 2; % sucsessfull burn
    elseif percent_consumption > 0.30
        burn_type_o = 1; % marginal burn
    else
        burn_type_o = 0; % unsucsessfull burn
    end

    % Determine EXP Burn Type
    if burned == 1
        burn_type_e = 2; % sucsessfull burn
    elseif pct_burn > 0
        burn_type_e = 1; % marginal burn
    else
        burn_type_e = 0; % no burn
    end

    %Im not sure why, but the exell spreadsheet has these burns like this
    if burn_no == 53
        burn_type_e = 0;
    elseif burn_no == 103 || burn_no ==104
        burn_type_e = 1;
    end

% seperates inaccurate and accurate FDS predictions
match_d = [burn_no,burn_type_e,burn_type_o,percent_consumption, thisline(5).*2.5,thisline(4),eff_spacing,eff_depth,mstcnt,pct_burn];
    if burn_type_e ~= burn_type_o
        no_match_data = vertcat(no_match_data,match_d);
    else
        match_data = vertcat(match_data,match_d);
    end
all_percent = vertcat(all_percent,match_d);
end

plotx = 1; % switch the x axis value
switch plotx
    case 1
        x_bound = [28,140];
        x_label = 'Burn Number';
        ld_nm = '';
    case 6
        x_bound = [-3,42];
        x_label = 'Slope (degrees)';
        ld_nm = '_slope';
    case 7
        x_bound = [15,45];
        x_label = 'Effective Horizontal Spacing (cm)';
        ld_nm = '_horzSp';
    case 8
        x_bound = [20,160];
        x_label = 'Effective Vertical Spacing (cm)';
        ld_nm = '_vertSp';
    case 9
        x_bound = [0.06,0.21];
        x_label = 'Moisture content';
        ld_nm = '_mc';
    case 10
        x_bound = [0,100];
        x_label = 'Burn % EXP';
        ld_nm = '_percentEXP';
end

plot_style
Title_Font_Size = 12;
Label_Font_Size = 12;
figure

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% plots fds burn type thresholds as lines and accurate fds runs as gray circles
p = plot([-5,200],[.60,.60],'-.k',[-5,200],[.30,.30],'--k',match_data(:,plotx),match_data(:,4),'o');
p(3).Color=[0.5 0.5 0.5];

hold on

% plots inaccurate FDS runs grouped according to their exp type
h = gscatter(no_match_data(:,plotx),no_match_data(:,4),no_match_data(:,2),"r","^",5);
h(1).Marker = 'diamond';
h(1).Color = 'g';
h(2).Marker = 'square';
h(2).Color = 'b';

xlim(x_bound) 
ylim([0,1])

title('Deep Fuel Beds FDS Fuel Consumption Fraction');
plotname ='USFS_Deep_Fuel_Beds_Fire_Spread_Threshold';
Ind_Title = x_label;
Dep_Title ='Fuel Energy Burned FDS (%)';

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend('good/marg fds bound','marg/bad fds bound', 'correct burns', ...
    'no burn exp','marginal burn exp','good burn exp','location','eastoutside');
set(lh,'FontName',Font_Name,'FontSize',8)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plotname,ld_nm]);

return