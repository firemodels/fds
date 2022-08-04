
%6/24/2022 Noelle Crump

close all
clear all

    % Runs from Matlab folder (Validation plots implementation)
params_loc = '../../Validation/Theobald_Hose_Stream/FDS_Input_Files/Build_Input_Files/';
outdir = '../../../out/Theobald_Hose_Stream/';
expdir = '../../../exp/Theobald_Hose_Stream/';
plot_dir = [pwd, '/../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Theobald_Hose_Stream/'];

%     % Runs from scripts folder (testing implementation)
% params_loc = '../../../Validation/Theobald_Hose_Stream/FDS_Input_Files/Build_Input_Files/';
% outdir = '../../../../out/Theobald_Hose_Stream/';
% expdir = '../../../../exp/Theobald_Hose_Stream/';
% plot_dir = [pwd, '/../../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Theobald_Hose_Stream/'];

chid_loc = [params_loc,'paramfile.csv'];
exp_loc = [expdir,'theobald_effect_1981_fds.csv'];

% get list of FDS CHIDs
import_chid = importdata(chid_loc, ',',1);
import_chid = import_chid.textdata;
chid_list = import_chid(:,2);
chid_list(1) = [];
s_Theobald = size(chid_list);

% import exp results
Theobald_import_exp = table2array(readtable(exp_loc,'PreserveVariableNames',true));
max_range_exp = Theobald_import_exp(:,8);
max_height_exp = Theobald_import_exp(:,6);
max_height_dist_exp = Theobald_import_exp(:,7);

% import exp parameters, fix NAN nozzle values
exp_noz = Theobald_import_exp(:,2);
exp_noz(isnan(exp_noz)) = 9;
exp_dia = Theobald_import_exp(:,3); 
exp_bar = Theobald_import_exp(:,4);
exp_deg = Theobald_import_exp(:,5);

% define lists for iteration
max_range_out = [];
max_height_out = [];
max_height_dist_out = [];

% iterate through data and build Plotting table
for runno = 1:s_Theobald(1)
    linefile = [chid_list{runno},'_line.csv'];
    out_import = importdata([outdir,linefile], ',', 2);

    % import fds data
    AMPUAx = out_import.data(:,1);
    AMPUA = out_import.data(:,2);
    AMPUA = AMPUA(~isnan(AMPUA));
    ZMAXz = out_import.data(:,3);
    ZMAX = out_import.data(:,4);
    ZMAX = ZMAX(~isnan(ZMAX));
    ZMAX_Xz = out_import.data(:,5);
    ZMAX_X = out_import.data(:,6);
    
    s_AMPUA = size(AMPUA);
    s_ZMAX = size(ZMAX);

    % find FDS max range
    f = s_AMPUA(1);
    while AMPUA(f) < 1e-3  % threashold value, list read end to start
        f = f-1;
        max_range = AMPUAx(f);
    end

    % find FDS max height
    f = s_ZMAX(1);
    while ZMAX(f) < 0.5 % threashold value, list read end to start
        f = f-1;
        max_height = ZMAXz(f); 
    end

    % get FDS max height dist
    max_height_dist = ZMAX_X(f);

    % build fds results lists
    max_range_out = vertcat(max_range_out,max_range);
    max_height_out = vertcat(max_height_out,max_height);
    max_height_dist_out = vertcat(max_height_dist_out,max_height_dist);

end


% Scale Diameter for plots
exp_dia1 = exp_dia*12;

% Remove NAN rows from max height and max dist arrays, make tables
TheoArr1 = horzcat(max_range_exp,max_range_out,exp_noz,exp_dia1,exp_bar,exp_deg);
TheoArr2 = horzcat(max_height_exp,max_height_out,exp_noz,exp_dia1,exp_bar,exp_deg);
TheoArr3 = horzcat(max_height_dist_exp,max_height_dist_out,exp_noz,exp_dia1,exp_bar,exp_deg);
TheoArr2(any(isnan(TheoArr2), 2), :) = [];
TheoArr3(any(isnan(TheoArr3), 2), :) = [];
T1 = array2table(TheoArr1,'VariableNames',{'max range exp','max range out','nozzle','diameter','operating pressure','firing angle'});
T2 = array2table(TheoArr2,'VariableNames',{'max height exp','max height out','nozzle','diameter','operating pressure','firing angle'});
T3 = array2table(TheoArr3,'VariableNames',{'max height dist exp','max height dist out','nozzle','diameter','operating pressure','firing angle'});

% start setting plot style parameters
plot_style

Title_Font_Size = 14;
Label_Font_Size = 14;

Plot_Width = 4.75;
Plot_Height = 4.75;
Plot_X  = 1.00;
Plot_Y  = 0.75;

Paper_Height = 6.0;
Paper_Width  = 6.5;

% Change if the ColorBar should be Pressure or Firing Angle
ColVar = 1;
switch ColVar
    case 1 
        Color_Variable = 'operating pressure';
        Ticks = [2.1,2.8,4.1,4.8,6.2];
        TickLabels = {'2.1','2.8','4.1','4.8','6.2'};
        ColorBar_Label = 'Operating Pressure (Pa) 1e5';
        filename_add = '';
    case 2 
        Color_Variable = 'firing angle';
        Ticks = [20,25,30,35,40,45];
        TickLabels = {'20','25','30','35','40','45'};
        ColorBar_Label = 'Firing Angle (deg)';
        filename_add = '_fa';
end

% ----------------- Plot 1 Max Range -----------------
plot_filename = 'Theobald_Hose_Stream_Max_Range';
Ind_Title = 'Measured Range EXP (m)';
Dep_Title = 'Predicted Range FDS (m)';

figure();
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

h = gscatter(TheoArr1(:,1),TheoArr1(:,2),TheoArr1(:,3),'c','v',5);
h(2).Marker = 'diamond';
h(3).Marker = 'square';
h(4).Marker = '^';
h(2).Color = 'm';
h(3).Color = 'r';
h(4).Color = 'g';

hold on

scatter(T1,'max range exp','max range out','SizeVariable','diameter','ColorVariable',Color_Variable);
c = colorbar('Ticks',Ticks,'TickLabels',TickLabels);
c.Label.String = ColorBar_Label;

lim = [15,65];
plot(lim,lim,'k-');
title('Theobald Hose Stream Max Range');
xlim(lim)
ylim(lim)

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend('Nozzle 6','Nozzle 7','Nozzle 9 (Rouse)','Nozzle 10','Nozzle Diameter','','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plot_filename,filename_add]);

hold off

% ----------------- Plot 2 Max Height -----------------
plot_filename = 'Theobald_Hose_Stream_Max_Height';
Ind_Title = 'Measured Height EXP (m)';
Dep_Title = 'Predicted Height FDS (m)';

figure();
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

h = gscatter(TheoArr2(:,1),TheoArr2(:,2),TheoArr2(:,3),'m','d',5);
h(2).Marker = 'square';
h(2).Color = 'r';

hold on

scatter(T2,'max height exp','max height out','SizeVariable','diameter','ColorVariable',Color_Variable);
c = colorbar('Ticks',Ticks,'TickLabels',TickLabels);
c.Label.String = ColorBar_Label;

lim = [0,18];
plot(lim,lim,'k-');
title('Theobald Hose Stream Max Height');
xlim(lim)
ylim(lim)

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend('Nozzle 7','Nozzle 9 (Rouse)','Nozzle Diameter','','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plot_filename,filename_add]);
hold off

% ----------------- Plot 3 Max Height Dist -----------------
plot_filename ='Theobald_Hose_Stream_Max_Height_Distance';
Ind_Title ='Measured Distance EXP (m)';
Dep_Title ='Predicted Distance FDS (m)';

figure();
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

h=gscatter(TheoArr3(:,1),TheoArr3(:,2),TheoArr3(:,3),'m','d',5);
h(2).Marker='square';
h(2).Color='r';

hold on
scatter(T3,'max height dist exp','max height dist out','SizeVariable','diameter','ColorVariable',Color_Variable);
c = colorbar('Ticks',Ticks,'TickLabels',TickLabels);
c.Label.String = ColorBar_Label;

lim = [4,40];
plot(lim,lim,'k-');
title('Theobald Hose Stream Max Height Distance');
xlim(lim)
ylim(lim)

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend('Nozzle 7','Nozzle 9 (Rouse)','Nozzle Diameter','','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plot_filename,filename_add]);
return