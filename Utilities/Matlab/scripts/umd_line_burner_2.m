% McDermott
% 1-9-2017
% umd_line_burner_2.m
%
% Plot combustion efficiency for UMD Line Burner cases.

close all
clear all

CH4_HOC = 50010.3475; % kJ/kg-CH4 from .out file

Lf_pts = 50;
Lf_min = 0;
Lf_max = 1;
Lf_dz  = (Lf_max-Lf_min)/(Lf_pts-1);
Lf_z   = Lf_min:Lf_dz:Lf_max;
Lf_fac = 0.99;

plot_style

expdir = '../../../exp/Submodules/macfp-db/Extinction/UMD_Line_Burner/Experimental_Data/';
outdir = '../../../out/UMD_Line_Burner/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';

% experimental results
EXP = importdata([expdir,'CH4_A_Data.csv'],',',1);
XO2   = EXP.data(:,find(strcmp(EXP.colheaders,'XO2')));
q_R   = EXP.data(:,find(strcmp(EXP.colheaders,'q_R')));
Chi_R = EXP.data(:,find(strcmp(EXP.colheaders,'Chi_R')));
S_XO2 = EXP.data(:,find(strcmp(EXP.colheaders,'S_XO2'))); % uncertainty
eta   = EXP.data(:,find(strcmp(EXP.colheaders,'eta')));
S_eta = EXP.data(:,find(strcmp(EXP.colheaders,'S_eta')));

EXP_Lf = importdata([expdir,'CH4_A_Lf_Data.csv'],',',1);
XO2_Lf = EXP_Lf.data(:,find(strcmp(EXP_Lf.colheaders,'XO2_Lf')));
Lf     = EXP_Lf.data(:,find(strcmp(EXP_Lf.colheaders,'Lf')));
S_Lf   = EXP_Lf.data(:,find(strcmp(EXP_Lf.colheaders,'S_Lf')));

% dx=1.25 cm results
HRR1 = importdata([outdir,'methane_XO2_ramp_dx_1p25cm_hrr.csv'],',',2);
DEV1 = importdata([outdir,'methane_XO2_ramp_dx_1p25cm_devc.csv'],',',2);

Time_FDS_1 = DEV1.data(:,find(strcmp(DEV1.colheaders,'Time')));
XO2_FDS_1 = DEV1.data(:,find(strcmp(DEV1.colheaders,'"XO2"')));
q_R_FDS_1 = 0.5*( DEV1.data(:,find(strcmp(DEV1.colheaders,'"qrad1"'))) + DEV1.data(:,find(strcmp(DEV1.colheaders,'"qrad2"'))) );
eta_FDS_1 = HRR1.data(:,find(strcmp(HRR1.colheaders,'HRR')));
MLR_FDS_1 = HRR1.data(:,find(strcmp(HRR1.colheaders,'MLR_FUEL')));
eta_FDS_1 = eta_FDS_1./(MLR_FDS_1*CH4_HOC);
CHI_R_FDS_1 = DEV1.data(:,find(strcmp(DEV1.colheaders,'"CHI_R"')));

% flame height
colLf01 = find(strcmp(DEV1.colheaders,'"Lf-01"'));
colLf50 = find(strcmp(DEV1.colheaders,'"Lf-50"'));
Lf_FDS_1 = zeros(1,length(Time_FDS_1));
for n=1:length(Time_FDS_1)
    hrrpul = DEV1.data(n,colLf01:colLf50);
    q_total = sum(hrrpul);
    q_sum = 0.;
    for i=1:Lf_pts
        q_sum = q_sum+hrrpul(i);
        if q_sum>Lf_fac*q_total
            Lf_FDS_1(n) = Lf_z(i);
            break;
        end
    end
end

% dx=0.625 cm results
HRR2 = importdata([outdir,'methane_XO2_ramp_dx_p625cm_hrr.csv'],',',2);
DEV2 = importdata([outdir,'methane_XO2_ramp_dx_p625cm_devc.csv'],',',2);

Time_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'Time')));
XO2_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'"XO2"')));
q_R_FDS_2 = 0.5*( DEV2.data(:,find(strcmp(DEV2.colheaders,'"qrad1"'))) + DEV2.data(:,find(strcmp(DEV2.colheaders,'"qrad2"'))) );
eta_FDS_2 = HRR2.data(:,find(strcmp(HRR2.colheaders,'HRR')));
MLR_FDS_2 = HRR2.data(:,find(strcmp(HRR2.colheaders,'MLR_FUEL')));
eta_FDS_2 = eta_FDS_2./(MLR_FDS_2*CH4_HOC);
CHI_R_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'"CHI_R"')));

% flame height
colLf01 = find(strcmp(DEV2.colheaders,'"Lf-01"'));
colLf50 = find(strcmp(DEV2.colheaders,'"Lf-50"'));
Lf_FDS_2 = zeros(1,length(Time_FDS_2));
for n=1:length(Time_FDS_2)
    hrrpul = DEV2.data(n,colLf01:colLf50);
    q_total = sum(hrrpul);
    q_sum = 0.;
    for i=1:Lf_pts
        q_sum = q_sum+hrrpul(i);
        if q_sum>Lf_fac*q_total
            Lf_FDS_2(n) = Lf_z(i);
            break;
        end
    end
end

% dx=0.3125 cm results
HRR3 = importdata([outdir,'methane_XO2_ramp_dx_p3125cm_hrr.csv'],',',2);
DEV3 = importdata([outdir,'methane_XO2_ramp_dx_p3125cm_devc.csv'],',',2);

Time_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'Time')));
XO2_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'"XO2"')));
q_R_FDS_3 = 0.5*( DEV3.data(:,find(strcmp(DEV3.colheaders,'"qrad1"'))) + DEV3.data(:,find(strcmp(DEV3.colheaders,'"qrad2"'))) );
eta_FDS_3 = HRR3.data(:,find(strcmp(HRR3.colheaders,'HRR')));
MLR_FDS_3 = HRR3.data(:,find(strcmp(HRR3.colheaders,'MLR_FUEL')));
eta_FDS_3 = eta_FDS_3./(MLR_FDS_3*CH4_HOC);
CHI_R_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'"CHI_R"')));

% flame height
colLf01 = find(strcmp(DEV3.colheaders,'"Lf-01"'));
colLf50 = find(strcmp(DEV3.colheaders,'"Lf-50"'));
Lf_FDS_3 = zeros(1,length(Time_FDS_3));
for n=1:length(Time_FDS_3)
    hrrpul = DEV3.data(n,colLf01:colLf50);
    q_total = sum(hrrpul);
    q_sum = 0.;
    for i=1:Lf_pts
        q_sum = q_sum+hrrpul(i);
        if q_sum>Lf_fac*q_total
            Lf_FDS_3(n) = Lf_z(i);
            break;
        end
    end
end

% compute CHI_R ramp
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

X_O2_max = 2.077972E-01;
X_O2_min = 0.1;
ramp_time = @(x) (x-X_O2_max)*60./(X_O2_min-X_O2_max); % this comes from linear XO2 ramp (see input file)
H(1)=plot(XO2,Chi_R); hold on
pw_ramp_time = [0,37,45,60];
pw_ramp_XO2  = [0.21,0.142,0.123,0];
pw_ramp_chir = [.245,.125,0,0];
H(2)=plot(pw_ramp_XO2,pw_ramp_chir);
H(3)=plot(XO2_FDS_1,CHI_R_FDS_1,'ksq');
H(4)=plot(XO2_FDS_2,CHI_R_FDS_2,'b^');
H(5)=plot(XO2_FDS_3,CHI_R_FDS_3,'m>');
axis([0.12 0.22 0 0.30 ])
xlabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('\chi_R','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(5,0.27,'Radiative Fraction Ramp','FontName',Font_Name,'FontSize',Title_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
lh=legend(H,'Exp','Linear Ramp','FDS {\itW/\deltax}=4','FDS {\itW/\deltax}=8','FDS {\itW/\deltax}=16','Location','NorthWest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
git_file=[outdir,'methane_XO2_ramp_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_Chi_r_ramp_check']);

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
H(2)=plot(XO2_FDS_1,eta_FDS_1,'-.');
H(3)=plot(XO2_FDS_2,eta_FDS_2,'--');
H(4)=plot(XO2_FDS_3,eta_FDS_3,'-');
axis([0.09 0.21 0 1.2 ])
xlabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('\eta','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
text(0.095,1.12,'Methane Combustion Efficiency','FontName',Font_Name,'FontSize',Title_Font_Size)

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
XO2_ramp = [X_O2_max,X_O2_min];
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(Time_ramp,XO2_ramp,'--o','MarkerSize',10); hold on
H(2)=plot(Time_FDS_1,XO2_FDS_1,'-*'); hold on
axis([0 60 0.1 0.21])
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
lh=legend(H,'Input Ramp','FDS','Location','NorthEast');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

git_file=[outdir,'methane_XO2_ramp_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_N2_ramp_check']);

% plot heat flux
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(XO2,q_R,'o'); hold on
H(2)=plot(XO2_FDS_1,q_R_FDS_1,'-.','LineWidth',3);
H(3)=plot(XO2_FDS_2,q_R_FDS_2,'--','LineWidth',3);
H(4)=plot(XO2_FDS_3,q_R_FDS_3,'-' ,'LineWidth',3);
axis([0.12 0.21 0 1.4])
xlabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\itq"}_R (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
lh=legend(H,'Exp','FDS {\itW/\deltax}=4','FDS {\itW/\deltax}=8','FDS {\itW/\deltax}=16','Location','NorthWest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

git_file=[outdir,'methane_XO2_ramp_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_rad_heat_flux']);

% plot flame height
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(XO2_Lf,Lf,'o'); hold on
subr=1:2:length(Lf);
h=errorbar(XO2_Lf(subr),Lf(subr),-S_Lf(subr),S_Lf(subr),'o'); hold on
set(h,'Color',steel_blue)

H(2)=plot(XO2_FDS_1,Lf_FDS_1,'-.');
H(3)=plot(XO2_FDS_2,Lf_FDS_2,'--');
H(4)=plot(XO2_FDS_3,Lf_FDS_3,'-');
axis([0.14 0.22 0.4 0.75])
xlabel('{\itX}_{O2}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\itL}_f (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
lh=legend(H,'Exp','FDS {\itW/\deltax}=4','FDS {\itW/\deltax}=8','FDS {\itW/\deltax}=16','Location','NorthWest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

git_file=[outdir,'methane_XO2_ramp_dx_1p25cm_git.txt'];
addverstr(gca,git_file,'linear');

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,'methane_flame_height']);


