% McDermott
% 1-9-2017
% umd_line_burner_2.m
%
% Plot combustion efficiency for UMD Line Burner cases.

close all
clear all

CH4_HOC = 50010.3475; % kJ/kg-CH4 from .out file

plot_style

expdir = '../../../exp/Submodules/macfp-db/Extinction/UMD_Line_Burner/Experimental_Data/';
outdir = '../../../out/UMD_Line_Burner/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';

% experimental results
EXP = importdata([expdir,'CH4_A_Data.csv'],',',1);
XO2   = EXP.data(:,find(strcmp(EXP.colheaders,'XO2')));
S_XO2 = EXP.data(:,find(strcmp(EXP.colheaders,'S_XO2'))); % uncertainty
eta   = EXP.data(:,find(strcmp(EXP.colheaders,'eta')));
S_eta = EXP.data(:,find(strcmp(EXP.colheaders,'S_eta')));
Chi_R = EXP.data(:,find(strcmp(EXP.colheaders,'Chi_R')));

% compute CHI_R ramp
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

ramp_time = @(x) (x-0.21)*60./(0.1-0.21); % this comes from linear XO2 ramp (see input file)
H(1)=plot(ramp_time(XO2),Chi_R); hold on
pw_ramp_time = [0,37,45,60];
pw_ramp_chir = [.245,.125,.03,0];
H(2)=plot(pw_ramp_time,pw_ramp_chir);
axis([0 60 0 0.30 ])
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('\chi_r','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(5,0.27,'Radiative Fraction Ramp','FontName',Font_Name,'FontSize',Title_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
lh=legend(H,'Exp','Linear Ramp','Location','NorthEast');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
git_file=[outdir,'methane_XO2_ramp_dx_p625cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_Chi_r_ramp_check']);

% dx=0.625 cm results
HRR = importdata([outdir,'methane_XO2_ramp_dx_p625cm_hrr.csv'],',',2);
DEV = importdata([outdir,'methane_XO2_ramp_dx_p625cm_devc.csv'],',',2);

Time_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'Time')));
XO2_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'"XO2"')));
eta_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'HRR')));
MLR_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'MLR_FUEL')));
eta_FDS = eta_FDS./(MLR_FDS*CH4_HOC);

% dx=1.25 cm results
HRR2 = importdata([outdir,'methane_XO2_ramp_dx_1p25cm_hrr.csv'],',',2);
DEV2 = importdata([outdir,'methane_XO2_ramp_dx_1p25cm_devc.csv'],',',2);

XO2_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'"XO2"')));
eta_FDS_2 = HRR2.data(:,find(strcmp(HRR2.colheaders,'HRR')));
MLR_FDS_2 = HRR2.data(:,find(strcmp(HRR2.colheaders,'MLR_FUEL')));
eta_FDS_2 = eta_FDS_2./(MLR_FDS_2*CH4_HOC);

% dx=0.3125 cm results
HRR3 = importdata([outdir,'methane_XO2_ramp_dx_p3125cm_hrr.csv'],',',2);
DEV3 = importdata([outdir,'methane_XO2_ramp_dx_p3125cm_devc.csv'],',',2);

XO2_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'"XO2"')));
eta_FDS_3 = HRR3.data(:,find(strcmp(HRR3.colheaders,'HRR')));
MLR_FDS_3 = HRR3.data(:,find(strcmp(HRR3.colheaders,'MLR_FUEL')));
eta_FDS_3 = eta_FDS_3./(MLR_FDS_3*CH4_HOC);

figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

steel_blue = [0 0.447 0.741];
H(1)=plot(XO2,eta,'.','MarkerSize',10); hold on
set(H(1),'Color',steel_blue)
subr= 1:10:length(XO2);
h=errorbar(XO2(subr),eta(subr),-S_eta(subr),S_eta(subr),-S_XO2(subr),S_XO2(subr),'.','MarkerSize',10); hold on
set(h,'Color',steel_blue)
H(2)=plot(XO2_FDS_2,eta_FDS_2,'k-.');
H(3)=plot(XO2_FDS,eta_FDS,'k--');
H(4)=plot(XO2_FDS_3,eta_FDS_3,'k-','LineWidth',2);
axis([0.09 0.21 0 1.2 ])
xlabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('\eta','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(0.095,1.1,'Methane Combustion Efficiency','FontName',Font_Name,'FontSize',Title_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

lh=legend(H,'Exp','FDS {\itW/\deltax}=4','FDS {\itW/\deltax}=8','FDS {\itW/\deltax}=16','Location','SouthEast');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

git_file=[outdir,'methane_XO2_ramp_dx_p3125cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_eta']);

% double check N2 ramp

% intended ramp:
Time_ramp = [0,60];
XO2_ramp = [0.21,0.10];
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(Time_ramp,XO2_ramp,'--o','MarkerSize',10); hold on
H(2)=plot(Time_FDS,XO2_FDS,'-'); hold on
axis([0 60 0.1 0.21])
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
lh=legend(H,'Input Ramp','FDS','Location','NorthEast');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

git_file=[outdir,'methane_XO2_ramp_dx_p625cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_N2_ramp_check']);







