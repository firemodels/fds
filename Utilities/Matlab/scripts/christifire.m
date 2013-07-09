% Matala
% 07-05-2013
% christifire.m
%
% This reads in and plots the data for CHRISTIFIRE test case.

close all
clear all

% Set global reaction rate parameters

addpath('../../Validation/CHRISTIFIRE/Experimental_Data')
addpath('../../Validation/CHRISTIFIRE/FDS_Output_Files')

close all
    
plot_style
set(gcf,'DefaultLineLineWidth',Line_Width)
set(gca,'FontName',Font_Name)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

skip_case = 0;
if ~exist('CHRISTIFIRE_S701_tga_N2_v1_devc.csv')
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_S701_tga_N2_v2_devc.csv')
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_S701_tga_air_v1_devc.csv')
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_S701_tga_air_v2_devc.csv')
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_S701_mcc_v1_devc.csv')
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_S701_mcc_v2_devc.csv')
    display('Error: File CHRISTIFIRE_S701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_I701_tga_N2_v1_devc.csv')
    display('Error: File CHRISTIFIRE_I701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_I701_tga_N2_v2_devc.csv')
    display('Error: File CHRISTIFIRE_I701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_I701_mcc_v1_devc.csv')
    display('Error: File CHRISTIFIRE_I701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_I701_mcc_v2_devc.csv')
    display('Error: File CHRISTIFIRE_I701_tga_N2_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_C701_cone_25_v1_devc.csv')
    display('Error: File CHRISTIFIRE_C701_cone_25_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_C701_cone_50_v1_devc.csv')
    display('Error: File CHRISTIFIRE_C701_cone_50_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_C701_cone_75_v1_devc.csv')
    display('Error: File CHRISTIFIRE_C701_cone_75_v1_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_C701_cone_25_v2_devc.csv')
    display('Error: File CHRISTIFIRE_C701_cone_25_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_C701_cone_50_v2_devc.csv')
    display('Error: File CHRISTIFIRE_C701_cone_50_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end
if ~exist('CHRISTIFIRE_C701_cone_75_v2_devc.csv')
    display('Error: File CHRISTIFIRE_C701_cone_75_v2_devc.csv does not exist. Skipping case.')
    skip_case = 1;
end

if skip_case
    return
end

% TGA (Cable sheath)---------------------------------------------
% read experimental results
[S_tga_N2] = csvread('CHRISTIFIRE_S701_tga_N2.csv',2);
[S_tga_air] = csvread('CHRISTIFIRE_S701_tga_air.csv',2);
% read FDS results
[S_tga_N2_v1] = csvread('CHRISTIFIRE_S701_tga_N2_v1_devc.csv',2);
[S_tga_air_v1] = csvread('CHRISTIFIRE_S701_tga_air_v1_devc.csv',2);
[S_tga_N2_v2] = csvread('CHRISTIFIRE_S701_tga_N2_v2_devc.csv',2);
[S_tga_air_v2] = csvread('CHRISTIFIRE_S701_tga_air_v2_devc.csv',2);

hf = figure(1);
hX = plot(S_tga_N2(:,1),S_tga_N2(:,2),'k');
hold on
plot(S_tga_N2_v1(:,2),S_tga_N2_v1(:,4)./max(S_tga_N2_v1(:,4)).*100, ...
    S_tga_N2_v2(:,2),S_tga_N2_v2(:,4)./max(S_tga_N2_v2(:,4)).*100)
plot(S_tga_air(:,1),S_tga_air(:,2),'k--')
plot(S_tga_N2_v1(:,2),S_tga_N2_v1(:,4)./max(S_tga_N2_v1(:,4)).*100, '--', ...
    S_tga_N2_v2(:,2),S_tga_N2_v2(:,4)./max(S_tga_N2_v2(:,4)).*100, '--')
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Temperature (\circC)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('Mass (%)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp (N2)','FDS v1 (N2)', 'FDS v2 (N2)', 'Exp (air)','FDS v1 (air)', 'FDS v2 (air)')
xlim([0,800])
ylim([0,100])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_S701_tga_N2_v1_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_S701_tga';
print(gcf,'-dpdf',plotname);

% MCC(Cable sheath)-------------------------------------------------
% read experimental results
[S_mcc] = csvread('CHRISTIFIRE_S701_mcc.csv',2);
% read FDS results
[S_mcc_v1] = csvread('CHRISTIFIRE_S701_mcc_v1_devc.csv',2);
[S_mcc_v2] = csvread('CHRISTIFIRE_S701_mcc_v2_devc.csv',2);

hf = figure(2);
hX = plot(S_mcc(:,1),S_mcc(:,2),'k');
hold on
plot(S_mcc_v1(:,2),S_mcc_v1(:,5), ...
    S_mcc_v2(:,2),S_mcc_v2(:,5))

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Temperature (\circC)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('HRR (kW/kg)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2')
xlim([100,600])
ylim([0,200])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_S701_mcc_v1_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_S701_mcc';
print(gcf,'-dpdf',plotname);


% TGA (Cable insulation)---------------------------------------------
% read experimental results
[I_tga_N2] = csvread('CHRISTIFIRE_I701_tga_N2.csv',2);
% read FDS results
[I_tga_N2_v1] = csvread('CHRISTIFIRE_I701_tga_N2_v1_devc.csv',2);
[I_tga_N2_v2] = csvread('CHRISTIFIRE_I701_tga_N2_v2_devc.csv',2);

hf = figure(3);
hX = plot(I_tga_N2(:,1),I_tga_N2(:,2),'k');
hold on
plot(I_tga_N2_v1(:,2),I_tga_N2_v1(:,4)./max(I_tga_N2_v1(:,4)).*100, ...
    I_tga_N2_v2(:,2),I_tga_N2_v2(:,4)./max(I_tga_N2_v2(:,4)).*100)

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Temperature (\circC)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('Mass (%)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2')
xlim([0,800])
ylim([0,100])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_I701_tga_N2_v1_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_I701_tga';
print(gcf,'-dpdf',plotname);

% MCC(Cable insulation)--------------------------------------------------
% read experimental results
[I_mcc] = csvread('CHRISTIFIRE_I701_mcc.csv',2);
% read FDS results
[I_mcc_v1] = csvread('CHRISTIFIRE_I701_mcc_v1_devc.csv',2);
[I_mcc_v2] = csvread('CHRISTIFIRE_I701_mcc_v2_devc.csv',2);

hf = figure(4);
hX = plot(I_mcc(:,1),I_mcc(:,2),'k');
hold on
plot(I_mcc_v1(:,2),I_mcc_v1(:,5), ...
    I_mcc_v2(:,2),I_mcc_v2(:,5))

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Temperature (\circC)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('HRR (kW/kg)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2')
xlim([100,600])
ylim([0,200])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_I701_mcc_v1_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_I701_mcc';
print(gcf,'-dpdf',plotname);

% Cone calorimeter results---------------------------------------------
% read experimental results
[C_cone_25] = csvread('CHRISTIFIRE_C701_cone_25.csv',2);
[C_cone_50] = csvread('CHRISTIFIRE_C701_cone_50.csv',2);
[C_cone_75] = csvread('CHRISTIFIRE_C701_cone_75.csv',2);
% read FDS results
[C_cone_25_v1] = csvread('CHRISTIFIRE_C701_cone_25_v1_devc.csv',2);
[C_cone_50_v1] = csvread('CHRISTIFIRE_C701_cone_50_v1_devc.csv',2);
[C_cone_75_v1] = csvread('CHRISTIFIRE_C701_cone_75_v1_devc.csv',2);
[C_cone_25_v2] = csvread('CHRISTIFIRE_C701_cone_25_v2_devc.csv',2);
[C_cone_50_v2] = csvread('CHRISTIFIRE_C701_cone_50_v2_devc.csv',2);
[C_cone_75_v2] = csvread('CHRISTIFIRE_C701_cone_75_v2_devc.csv',2);

%plot HRR 50 kW/m2
hf = figure(5);
hX = plot(C_cone_50(:,1),C_cone_50(:,2),'k');
hold on
plot(C_cone_50_v1(:,1),C_cone_50_v1(:,3)./0.01, ...
    C_cone_50_v2(:,1),C_cone_50_v2(:,3)./0.01) %scaled with the area

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Time (s)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('HRR (kW/m^2)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2')
xlim([0,1500])
ylim([0,350])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_cone_50_v1_devc_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_hrr_50';
print(gcf,'-dpdf',plotname);

%plot HRR 25 and 75 kW/m2
hf = figure(6);
hX = plot(C_cone_25(:,1),C_cone_25(:,2),'k');
hold on
plot(C_cone_25_v1(:,1),C_cone_25_v1(:,3)./0.01, ...
    C_cone_25_v2(:,1),C_cone_25_v2(:,3)./0.01) %scaled with the area
plot(C_cone_75(:,1),C_cone_75(:,2),'k--')
plot(C_cone_75_v1(:,1),C_cone_75_v1(:,3)./0.01,'--', ...
    C_cone_75_v2(:,1),C_cone_75_v2(:,3)./0.01,'--') %scaled with the area

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Time (s)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('HRR (kW/m^2)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp 25','FDS v1 25', 'FDS v2 25',...
    'Exp 75','FDS v1 75', 'FDS v2 75')
xlim([0,2500])
ylim([0,400])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_cone_50_v1_devc_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_hrr_25_75';
print(gcf,'-dpdf',plotname);

%plot MLR 50 kW/m2
hf = figure(7);
hX = plot(C_cone_50(:,1),C_cone_50(:,3),'k');
hold on
plot(C_cone_50_v1(:,1),-gradient(C_cone_50_v1(:,2),C_cone_50_v1(:,1)), ...
    C_cone_50_v2(:,1),-gradient(C_cone_50_v2(:,2),C_cone_50_v2(:,1))) %scaled with the area

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Time (s)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('MLR (kg/sm^2)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2')
xlim([0,1500])
ylim([0,0.016])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_cone_50_v1_devc_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_mlr_50';
print(gcf,'-dpdf',plotname);

%plot MLR 25 and 75 kW/m2
hf = figure(8);
hX = plot(C_cone_25(:,1),C_cone_25(:,3),'k');
hold on
plot(C_cone_25_v1(:,1),-gradient(C_cone_25_v1(:,2),C_cone_25_v1(:,1)), ...
    C_cone_25_v2(:,1),-gradient(C_cone_25_v2(:,2),C_cone_25_v2(:,1))) %scaled with the area
plot(C_cone_75(:,1),C_cone_75(:,3),'k--')
plot(C_cone_75_v1(:,1),-gradient(C_cone_75_v1(:,2),C_cone_75_v1(:,1)),'--', ...
    C_cone_75_v2(:,1),-gradient(C_cone_75_v2(:,2),C_cone_75_v2(:,1)),'--') %scaled with the area

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Time (s)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('MLR (kg/sm^2)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp 25','FDS v1 25', 'FDS v2 25',...
    'Exp 75','FDS v1 75', 'FDS v2 75')
xlim([0,2500])
ylim([0,0.022])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_cone_50_v1_devc_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_mlr_25_75';
print(gcf,'-dpdf',plotname);

%plot EHC (effective heat of combustion) 50 kW/m2
hf = figure(9);
hX = plot(C_cone_50(:,1),C_cone_50(:,4),'k');
hold on
plot(C_cone_50_v1(:,1),C_cone_50_v1(:,3)./0.01./-gradient(C_cone_50_v1(:,2),C_cone_50_v1(:,1)), ...
    C_cone_50_v2(:,1),C_cone_50_v2(:,3)./0.01./-gradient(C_cone_50_v2(:,2),C_cone_50_v2(:,1))) %scaled with the area

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Time (s)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('EHC (kJ/kg)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp','FDS v1', 'FDS v2')
xlim([0,1500])
ylim([0,40000])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_cone_50_v1_devc_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_ehc_50';
print(gcf,'-dpdf',plotname);

%plot MLR 25 and 75 kW/m2
hf = figure(8);
hX = plot(C_cone_25(:,1),C_cone_25(:,4),'k');
hold on
plot(C_cone_25_v1(:,1),C_cone_25_v1(:,3)./0.01./-gradient(C_cone_25_v1(:,2),C_cone_25_v1(:,1)), ...
    C_cone_25_v2(:,1),C_cone_25_v2(:,3)./0.01./-gradient(C_cone_25_v2(:,2),C_cone_25_v2(:,1))) %scaled with the area
plot(C_cone_75(:,1),C_cone_75(:,4),'k--')
plot(C_cone_75_v1(:,1),C_cone_75_v1(:,3)./0.01./-gradient(C_cone_75_v1(:,2),C_cone_75_v1(:,1)), '--', ...
    C_cone_75_v2(:,1),C_cone_75_v2(:,3)./0.01./-gradient(C_cone_75_v2(:,2),C_cone_75_v2(:,1)),'--') %scaled with the area

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Scat_Plot_X,Scat_Plot_Y,Scat_Plot_Width,Scat_Plot_Height])
set(hf(1),'DefaultLineLineWidth',Line_Width)
xlabel('Time (s)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
ylabel('EHC (kJ/kg)','Interpreter','tex','FontSize',Scat_Label_Font_Size)
legend('Exp 25','FDS v1 25', 'FDS v2 25',...
    'Exp 75','FDS v1 75', 'FDS v2 75')
xlim([0,2500])
ylim([0,40000])

% add SVN if file is available

svn_file = '../../Validation/CHRISTIFIRE/FDS_Output_Files/CHRISTIFIRE_cone_50_v1_devc_svn.txt';

if exist(svn_file,'file')
   SVN = importdata(svn_file);
   x_lim = get(gca,'XLim');
   y_lim = get(gca,'YLim');
   X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
   Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
   text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
       'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Scat_Paper_Width Scat_Paper_Height]);
set(gcf,'PaperPosition',[0 0 Scat_Paper_Width Scat_Paper_Height]);
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/CHRISTIFIRE/CHRISTIFIRE_C701_ehc_25_75';
print(gcf,'-dpdf',plotname);

clear all
close all