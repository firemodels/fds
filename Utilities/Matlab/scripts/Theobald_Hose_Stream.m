%6/24/2022 Noelle Crump

close all
clear all

params_loc = '../../Validation/Theobald_Hose_Stream/FDS_Input_Files/Build_Input_Files/';
outdir = '../../../out/Theobald_Hose_Stream/';
expdir = '../../../exp/Theobald_Hose_Stream/';
plot_dir = ['../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Theobald_Hose_Stream/'];

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
    AMPUAx = out_import.data(:,find(strcmp(out_import.colheaders,'AMPUA-x')));
    AMPUA = out_import.data(:,find(strcmp(out_import.colheaders,'AMPUA')));
    AMPUA = AMPUA(~isnan(AMPUA));
    ZMAXz = out_import.data(:,find(strcmp(out_import.colheaders,'ZMAX-z')));
    ZMAX = out_import.data(:,find(strcmp(out_import.colheaders,'ZMAX')));
    ZMAX = ZMAX(~isnan(ZMAX));
    ZMAX_Xz = out_import.data(:,find(strcmp(out_import.colheaders,'ZMAX_X-z')));
    ZMAX_X = out_import.data(:,find(strcmp(out_import.colheaders,'ZMAX_X')));
    
    % find FDS max range
    for f=length(AMPUA):-1:1
        max_range = AMPUAx(f);
        if AMPUA(f) > 1
            break
        end
    end

    % find FDS max height
    for f=length(ZMAX):-1:1
        max_height = ZMAXz(f);
        if ZMAX(f) > 1
            break
        end
    end

    % get max height dist
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

% ----------------- Plot 1 Max Range -----------------
plot_filename = 'Theobald_Hose_Stream_Max_Range';
Ind_Title = 'Measured Max Range EXP (m)';
Dep_Title = 'Predicted Max Range FDS (m)';

figure();
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

range_6 = find(T1.('nozzle')==6);
h(1)=scatter(T1(range_6,:).('max range exp'),T1(range_6,:).('max range out'),'filled');
h(1).SizeData=T1{range_6,'diameter'}*.2;
h(1).CData=[0.8500 0.3250 0.0980];

hold on

range_7 = find(T1.('nozzle')==7);
h(2)=scatter(T1(range_7,:).('max range exp'),T1(range_7,:).('max range out'),'filled');
h(2).SizeData=T1{range_7,'diameter'}*.2;
h(2).CData=[0 0.4470 0.7410];

range_9 = find(T1.('nozzle')==9);
h(3)=scatter(T1(range_9,:).('max range exp'),T1(range_9,:).('max range out'),'filled');
h(3).SizeData=T1{range_9,'diameter'}*.2;
h(3).CData=[0.9290 0.6940 0.1250];

range_10 = find(T1.('nozzle')==10);
h(4)=scatter(T1(range_10,:).('max range exp'),T1(range_10,:).('max range out'),'filled');
h(4).SizeData=T1{range_10,'diameter'}*.2;
h(4).CData=[0.4940 0.1840 0.5560];

set(gca,'Box','on')

lim = [15,65];
plot(lim,lim,'k-');
% title('Theobald Hose Stream Max Range');
xlim(lim)
ylim(lim)

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend('Nozzle 6','Nozzle 7','Nozzle 9 (Rouse)','Nozzle 10','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git version if file is available
Git_Filename = [outdir,'Theobald_Test_0_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plot_filename]);

hold off

% ----------------- Plot 2 Max Height -----------------
plot_filename = 'Theobald_Hose_Stream_Max_Height';
Ind_Title = 'Measured Max Height EXP (m)';
Dep_Title = 'Predicted Max Height FDS (m)';

figure();
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

range_7 = find(T2.('nozzle')==7);
h(1)=scatter(T2(range_7,:).('max height exp'),T2(range_7,:).('max height out'),'filled');
h(1).SizeData=T2{range_7,'diameter'}*.2;
h(1).CData=[0 0.4470 0.7410];

hold on

range_9 = find(T2.('nozzle')==9);
h(2)=scatter(T2(range_9,:).('max height exp'),T2(range_9,:).('max height out'),'filled');
h(2).SizeData=T2{range_9,'diameter'}*.2;
h(2).CData=[0.9290 0.6940 0.1250];

set(gca,'Box','on')

hold on

lim = [0,18];
plot(lim,lim,'k-');
% title('Theobald Hose Stream Max Height');
xlim(lim)
ylim(lim)

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend('Nozzle 7','Nozzle 9 (Rouse)','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git version if file is available
Git_Filename = [outdir,'Theobald_Test_0_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plot_filename]);
hold off

% ----------------- Plot 3 Max Height Dist -----------------
plot_filename ='Theobald_Hose_Stream_Max_Height_Distance';
Ind_Title ='Measured Max Height Distance EXP (m)';
Dep_Title ='Predicted Max Height Distance FDS (m)';

figure();
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

range_7 = find(T3.('nozzle')==7);
h(1)=scatter(T3(range_7,:).('max height dist exp'),T3(range_7,:).('max height dist out'),'filled');
h(1).SizeData=T3{range_7,'diameter'}*.2;
h(1).CData=[0 0.4470 0.7410];

hold on

range_9 = find(T3.('nozzle')==9);
h(2)=scatter(T3(range_9,:).('max height dist exp'),T3(range_9,:).('max height dist out'),'filled');
h(2).SizeData=T3{range_9,'diameter'}*.2;
h(2).CData=[0.9290 0.6940 0.1250];

set(gca,'Box','on')

hold on

lim = [5,45];
plot(lim,lim,'k-');
% title('Theobald Hose Stream Max Height Distance');
xlim(lim)
ylim(lim)

xlabel(Ind_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel(Dep_Title,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh = legend(h(1:2),'Nozzle 7','Nozzle 9 (Rouse)','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add Git version if file is available
Git_Filename = [outdir,'Theobald_Test_0_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[plot_dir,plot_filename]);

