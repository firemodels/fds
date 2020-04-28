% Matala
% 07-05-2013
% christifire.m
%
% This reads in and plots the data for CHRISTIFIRE test case.

close all
clear all

disp('christifire ...')

% Set global reaction rate parameters

expdir = '../../../exp/CHRISTIFIRE/';
outdir = '../../../out/CHRISTIFIRE/';

skip_case = 0;
if ~exist([outdir,'CHRISTIFIRE_S701_tga_N2_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_S701_tga_N2_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_S701_tga_N2_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_S701_tga_air_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_S701_tga_air_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_S701_tga_air_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_S701_tga_air_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_S701_mcc_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_S701_mcc_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_S701_mcc_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_S701_mcc_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_I701_tga_N2_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_I701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_I701_tga_N2_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_I701_tga_N2_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_I701_mcc_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_I701_mcc_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_I701_mcc_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_I701_mcc_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_C701_cone_25_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_C701_cone_25_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_C701_cone_50_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_C701_cone_50_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_C701_cone_75_v1_devc.csv'])
    display('Error: File CHRISTIFIRE_C701_cone_75_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_C701_cone_25_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_C701_cone_25_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_C701_cone_50_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_C701_cone_50_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist([outdir,'CHRISTIFIRE_C701_cone_75_v2_devc.csv'])
    display('Error: File CHRISTIFIRE_C701_cone_75_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end

if skip_case
    return
end

% TGA (Cable sheath)---------------------------------------------
% read experimental results
[S_tga_N2] = csvread([expdir,'CHRISTIFIRE_S701_tga_N2.csv'],2);
[S_tga_air] = csvread([expdir,'CHRISTIFIRE_S701_tga_air.csv'],2);
% read FDS results
[S_tga_N2_v1] = csvread([outdir,'CHRISTIFIRE_S701_tga_N2_v1_devc.csv'],2);
[S_tga_air_v1] = csvread([outdir,'CHRISTIFIRE_S701_tga_air_v1_devc.csv'],2);
[S_tga_N2_v2] = csvread([outdir,'CHRISTIFIRE_S701_tga_N2_v2_devc.csv'],2);
[S_tga_air_v2] = csvread([outdir,'CHRISTIFIRE_S701_tga_air_v2_devc.csv'],2);

figure

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(S_tga_N2(:,1),S_tga_N2(:,2),'k');
hold on
plot(S_tga_N2_v1(:,2),S_tga_N2_v1(:,4)./max(S_tga_N2_v1(:,4)).*100, ...
    S_tga_N2_v2(:,2),S_tga_N2_v2(:,4)./max(S_tga_N2_v2(:,4)).*100)
plot(S_tga_air(:,1),S_tga_air(:,2),'k--')
plot(S_tga_air_v1(:,2),S_tga_air_v1(:,4)./max(S_tga_air_v1(:,4)).*100, '--', ...
    S_tga_air_v2(:,2),S_tga_air_v2(:,4)./max(S_tga_air_v2(:,4)).*100, '--')

xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Mass (%)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp (N2)','FDS v1 (N2)', 'FDS v2 (N2)', 'Exp (air)','FDS v1 (air)', 'FDS v2 (air)');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlim([0,800])
ylim([0,100])
Title_Scale_X = .05;
Title_Scale_Y = .05;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 sheath TGA','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_S701_tga_N2_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_S701_tga';
print(gcf,'-dpdf',plotname);

%return
close all

% MCC(Cable sheath)-------------------------------------------------
% read experimental results
[S_mcc] = csvread([expdir,'CHRISTIFIRE_S701_mcc.csv'],2);
% read FDS results
[S_mcc_v1] = csvread([outdir,'CHRISTIFIRE_S701_mcc_v1_devc.csv'],2);
[S_mcc_v2] = csvread([outdir,'CHRISTIFIRE_S701_mcc_v2_devc.csv'],2);

figure

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(S_mcc(:,1),S_mcc(:,2),'k');
hold on
plot(S_mcc_v1(:,2),S_mcc_v1(:,5), '--',  ...
    S_mcc_v2(:,2),S_mcc_v2(:,5), '--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('HRR (kW/kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([100,600])
ylim([0,250])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 sheath MCC','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_S701_mcc_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_S701_mcc';
print(gcf,'-dpdf',plotname);

%return
close all

% TGA (Cable insulation)---------------------------------------------
% read experimental results
[I_tga_N2] = csvread([expdir,'CHRISTIFIRE_I701_tga_N2.csv'],2);
% read FDS results
[I_tga_N2_v1] = csvread([outdir,'CHRISTIFIRE_I701_tga_N2_v1_devc.csv'],2);
[I_tga_N2_v2] = csvread([outdir,'CHRISTIFIRE_I701_tga_N2_v2_devc.csv'],2);

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(I_tga_N2(:,1),I_tga_N2(:,2),'k');
hold on
plot(I_tga_N2_v1(:,2),I_tga_N2_v1(:,4)./max(I_tga_N2_v1(:,4)).*100, '--', ...
    I_tga_N2_v2(:,2),I_tga_N2_v2(:,4)./max(I_tga_N2_v2(:,4)).*100, '--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Mass (%)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,800])
ylim([0,100])
Title_Scale_X = .05;
Title_Scale_Y = .05;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 insulation TGA','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_I701_tga_N2_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_I701_tga';
print(gcf,'-dpdf',plotname);

%return
close all

% MCC(Cable insulation)--------------------------------------------------
% read experimental results
[I_mcc] = csvread([expdir,'CHRISTIFIRE_I701_mcc.csv'],2);
% read FDS results
[I_mcc_v1] = csvread([outdir,'CHRISTIFIRE_I701_mcc_v1_devc.csv'],2);
[I_mcc_v2] = csvread([outdir,'CHRISTIFIRE_I701_mcc_v2_devc.csv'],2);

figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(I_mcc(:,1),I_mcc(:,2),'k');
hold on
plot(I_mcc_v1(:,2),I_mcc_v1(:,5), '--',  ...
    I_mcc_v2(:,2),I_mcc_v2(:,5), '--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Temperature (\circC)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('HRR (kW/kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([100,600])
ylim([0,250])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 insulation MCC','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_I701_mcc_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_I701_mcc';
print(gcf,'-dpdf',plotname);

%return
close all

% Cone calorimeter results---------------------------------------------
% read experimental results
[C_cone_25] = csvread([expdir,'CHRISTIFIRE_C701_cone_25.csv'],2);
[C_cone_50] = csvread([expdir,'CHRISTIFIRE_C701_cone_50.csv'],2);
[C_cone_75] = csvread([expdir,'CHRISTIFIRE_C701_cone_75.csv'],2);
% read FDS results
[C_cone_25_v1] = csvread([outdir,'CHRISTIFIRE_C701_cone_25_v1_devc.csv'],2);
[C_cone_50_v1] = csvread([outdir,'CHRISTIFIRE_C701_cone_50_v1_devc.csv'],2);
[C_cone_75_v1] = csvread([outdir,'CHRISTIFIRE_C701_cone_75_v1_devc.csv'],2);
[C_cone_25_v2] = csvread([outdir,'CHRISTIFIRE_C701_cone_25_v2_devc.csv'],2);
[C_cone_50_v2] = csvread([outdir,'CHRISTIFIRE_C701_cone_50_v2_devc.csv'],2);
[C_cone_75_v2] = csvread([outdir,'CHRISTIFIRE_C701_cone_75_v2_devc.csv'],2);

%plot HRR 50 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_50(:,1),C_cone_50(:,2),'k');
hold on
plot(C_cone_50_v1(:,1),C_cone_50_v1(:,3)./0.01, '--',  ...
    C_cone_50_v2(:,1),C_cone_50_v2(:,3)./0.01, '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('HRR (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,1500])
ylim([0,350])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 50','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_50_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_hrr_50';
print(gcf,'-dpdf',plotname);

%return
close all

%plot HRR 25 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_25(:,1),C_cone_25(:,2),'k');
hold on
plot(C_cone_25_v1(:,1),C_cone_25_v1(:,3)./0.01, '--', ...
    C_cone_25_v2(:,1),C_cone_25_v2(:,3)./0.01, '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('HRR (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2', 'Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,2500])
ylim([0,250])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 25','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_25_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_hrr_25';
print(gcf,'-dpdf',plotname);

%return
close all

%plot HRR 75 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_75(:,1),C_cone_75(:,2),'k');
hold on
plot(C_cone_75_v1(:,1),C_cone_75_v1(:,3)./0.01, '--', ...
    C_cone_75_v2(:,1),C_cone_75_v2(:,3)./0.01, '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('HRR (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,1500])
ylim([0,450])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 75','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_75_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_hrr_75';
print(gcf,'-dpdf',plotname);

%return
close all

%plot MLR 50 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_50(:,1),C_cone_50(:,3),'k');
hold on
plot(C_cone_50_v1(:,1),-gradient(C_cone_50_v1(:,2),C_cone_50_v1(:,1)), '--', ...
    C_cone_50_v2(:,1),-gradient(C_cone_50_v2(:,2),C_cone_50_v2(:,1)), '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('MLR (kg/sm^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,1500])
ylim([0,0.016])
Title_Scale_X = .05;
Title_Scale_Y = .05;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 50','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_50_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_mlr_50';
print(gcf,'-dpdf',plotname);

%return
close all

%plot MLR 25 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_25(:,1),C_cone_25(:,3),'k');
hold on
plot(C_cone_25_v1(:,1),-gradient(C_cone_25_v1(:,2),C_cone_25_v1(:,1)), '--', ...
    C_cone_25_v2(:,1),-gradient(C_cone_25_v2(:,2),C_cone_25_v2(:,1)), '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('MLR (kg/sm^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,2500])
ylim([0,0.015])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 25','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_25_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_mlr_25';
print(gcf,'-dpdf',plotname);

%return
close all

%plot MLR 75 kW/m2
figure(10)
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_75(:,1),C_cone_75(:,3),'k');
hold on
plot(C_cone_75_v1(:,1),-gradient(C_cone_75_v1(:,2),C_cone_75_v1(:,1)), '--', ...
    C_cone_75_v2(:,1),-gradient(C_cone_75_v2(:,2),C_cone_75_v2(:,1)), '--'); %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('MLR (kg/sm^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,1500])
ylim([0,0.03])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 75','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_75_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_mlr_75';
print(gcf,'-dpdf',plotname);

%return
close all

%plot EHC (effective heat of combustion) 50 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_50(:,1),C_cone_50(:,4),'k');
hold on
plot(C_cone_50_v1(:,1),C_cone_50_v1(:,3)./0.01./-gradient(C_cone_50_v1(:,2),C_cone_50_v1(:,1)), '--', ...
    C_cone_50_v2(:,1),C_cone_50_v2(:,3)./0.01./-gradient(C_cone_50_v2(:,2),C_cone_50_v2(:,1)), '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('EHC (kJ/kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2','Location', 'NorthWest')

xlim([0,1500])
ylim([0,40000])
Title_Scale_X = .05;
Title_Scale_Y = .05;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 50','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_50_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_ehc_50';
print(gcf,'-dpdf',plotname);

%return
close all

%plot EHC 25 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_25(:,1),C_cone_25(:,4),'k');
hold on
plot(C_cone_25_v1(:,1),C_cone_25_v1(:,3)./0.01./-gradient(C_cone_25_v1(:,2),C_cone_25_v1(:,1)), '--', ...
    C_cone_25_v2(:,1),C_cone_25_v2(:,3)./0.01./-gradient(C_cone_25_v2(:,2),C_cone_25_v2(:,1)), '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('EHC (kJ/kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,2500])
ylim([0,40000])
Title_Scale_X = .05;
Title_Scale_Y = .95;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 25','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_25_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_ehc_25';
print(gcf,'-dpdf',plotname);

%return
close all

%plot EHC 75 kW/m2
figure
plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(C_cone_75(:,1),C_cone_75(:,4),'k');
hold on
plot(C_cone_75_v1(:,1),C_cone_75_v1(:,3)./0.01./-gradient(C_cone_75_v1(:,2),C_cone_75_v1(:,1)), '--', ...
    C_cone_75_v2(:,1),C_cone_75_v2(:,3)./0.01./-gradient(C_cone_75_v2(:,2),C_cone_75_v2(:,1)), '--') %scaled with the area

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('EHC (kJ/kg)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
lh=legend('Exp','FDS v1', 'FDS v2','Location', 'East');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

xlim([0,1500])
ylim([0,40000])
Title_Scale_X = .05;
Title_Scale_Y = .05;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
X_Title_Position = x_lim(1)+Title_Scale_X*(x_lim(2)-x_lim(1));
Y_Title_Position = y_lim(1)+Title_Scale_Y*(y_lim(2)-y_lim(1));
text(X_Title_Position,Y_Title_Position,'CHRISTIFIRE cable 701 cone 75','FontName',Font_Name,'FontSize',Title_Font_Size)

% add VerStr if file is available
git_file = [outdir,'CHRISTIFIRE_C701_cone_75_v1_git.txt'];
addverstr(gca,git_file,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_ehc_75';
print(gcf,'-dpdf',plotname);

%return

clear all
close all
