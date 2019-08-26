% McDermott
% 1-9-2017
% umd_line_burner_2.m
%
% Plot combustion efficiency for UMD Line Burner cases.

close all
clear all

Lf_pts = 50;
Lf_min = 0;
Lf_max = 1;
Lf_dz  = (Lf_max-Lf_min)/(Lf_pts-1);
Lf_z   = Lf_min:Lf_dz:Lf_max;
Lf_fac = 0.99;
Lf_dt  = 10;

plot_style

expdir = '../../../exp/Submodules/macfp-db/Extinction/UMD_Line_Burner/Experimental_Data/';
outdir = '../../../out/UMD_Line_Burner/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';

exp_fname    = {'CH4_A_Data.csv','C3H8_A_Data.csv'};
exp_Lf_fname = {'CH4_A_Lf_Data.csv','C3H8_A_Lf_Data.csv'};
fuel_name    = {'methane','propane'};
Fuel_name    = {'Methane','Propane'};
fuel_hoc     = [50010.3475,46334.6246]; % from .out file
git_tag_ext  = '_2step_XO2_ramp_dx_p3125cm_git.txt';

line_fmt = {'bo-','ro-','ko-','b^-','r^-','k^-'};
key_fmt  = {'1step, {\it W/dx}=4','1step, {\it W/dx}=8','1step, {\it W/dx}=16','2step, {\it W/dx}=4','2step, {\it W/dx}=8','2step, {\it W/dx}=16'};

i_fuel = 2;

% experimental results
EXP = importdata([expdir,exp_fname{i_fuel}],',',1);
XO2   = EXP.data(:,find(strcmp(EXP.colheaders,'XO2')));
q_R   = EXP.data(:,find(strcmp(EXP.colheaders,'q_R')));
Chi_R = EXP.data(:,find(strcmp(EXP.colheaders,'Chi_R')));
S_XO2 = EXP.data(:,find(strcmp(EXP.colheaders,'S_XO2'))); % uncertainty
eta   = EXP.data(:,find(strcmp(EXP.colheaders,'eta')));
S_eta = EXP.data(:,find(strcmp(EXP.colheaders,'S_eta')));

EXP_Lf = importdata([expdir,exp_Lf_fname{i_fuel}],',',1);
XO2_Lf = EXP_Lf.data(:,find(strcmp(EXP_Lf.colheaders,'XO2_Lf')));
Lf     = EXP_Lf.data(:,find(strcmp(EXP_Lf.colheaders,'Lf')));
S_Lf   = EXP_Lf.data(:,find(strcmp(EXP_Lf.colheaders,'S_Lf')));

is_case_1_run=1; if ~exist([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_1p25cm_hrr.csv']); is_case_1_run=0; end
is_case_2_run=1; if ~exist([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_p625cm_hrr.csv']); is_case_2_run=0; end
is_case_3_run=1; if ~exist([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_p3125cm_hrr.csv']); is_case_3_run=0; end

is_case_4_run=1; if ~exist([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_1p25cm_hrr.csv']); is_case_4_run=0; end
is_case_5_run=1; if ~exist([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_p625cm_hrr.csv']); is_case_5_run=0; end
is_case_6_run=1; if ~exist([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_p3125cm_hrr.csv']); is_case_6_run=0; end

if is_case_1_run % case 1 if
    % dx=1.25 cm results
    HRR1 = importdata([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_1p25cm_hrr.csv'],',',2);
    DEV1 = importdata([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_1p25cm_devc.csv'],',',2);

    Time_FDS_1 = DEV1.data(:,find(strcmp(DEV1.colheaders,'Time')));
    XO2_FDS_1 = DEV1.data(:,find(strcmp(DEV1.colheaders,'"XO2"')));
    q_R_FDS_1 = 0.5*( DEV1.data(:,find(strcmp(DEV1.colheaders,'"qrad1"'))) + DEV1.data(:,find(strcmp(DEV1.colheaders,'"qrad2"'))) );
    HRR_FDS_1 = HRR1.data(:,find(strcmp(HRR1.colheaders,'HRR')));
    MLR_FDS_1 = HRR1.data(:,find(strcmp(HRR1.colheaders,'MLR_FUEL')));
    eta_FDS_1 = HRR_FDS_1./(MLR_FDS_1*fuel_hoc(i_fuel));
    CHI_R_FDS_1 = DEV1.data(:,find(strcmp(DEV1.colheaders,'"CHI_R"')));
    QR_FDS_1 = abs(HRR1.data(:,find(strcmp(HRR1.colheaders,'Q_RADI'))));
    GLOB_CHI_R_1 = QR_FDS_1./HRR_FDS_1;

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
    Lf_tmp = Lf_FDS_1;
    for n=1:length(Time_FDS_1)
        indx_range = [find(Time_FDS_1>(Time_FDS_1(n)-Lf_dt),1):n];
        Lf_FDS_1(n) = mean(Lf_tmp(indx_range));
    end
end % case 1 if

if is_case_2_run % case 2 if
    % dx=0.625 cm results
    HRR2 = importdata([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_p625cm_hrr.csv'],',',2);
    DEV2 = importdata([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_p625cm_devc.csv'],',',2);

    Time_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'Time')));
    XO2_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'"XO2"')));
    q_R_FDS_2 = 0.5*( DEV2.data(:,find(strcmp(DEV2.colheaders,'"qrad1"'))) + DEV2.data(:,find(strcmp(DEV2.colheaders,'"qrad2"'))) );
    HRR_FDS_2 = HRR2.data(:,find(strcmp(HRR2.colheaders,'HRR')));
    MLR_FDS_2 = HRR2.data(:,find(strcmp(HRR2.colheaders,'MLR_FUEL')));
    eta_FDS_2 = HRR_FDS_2./(MLR_FDS_2*fuel_hoc(i_fuel));
    CHI_R_FDS_2 = DEV2.data(:,find(strcmp(DEV2.colheaders,'"CHI_R"')));
    QR_FDS_2 = abs(HRR2.data(:,find(strcmp(HRR2.colheaders,'Q_RADI'))));
    GLOB_CHI_R_2 = QR_FDS_2./HRR_FDS_2;

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
    Lf_tmp = Lf_FDS_2;
    for n=1:length(Time_FDS_2)
        indx_range = [find(Time_FDS_2>(Time_FDS_2(n)-Lf_dt),1):n];
        Lf_FDS_2(n) = mean(Lf_tmp(indx_range));
    end
end % case 2 if

if is_case_3_run % case 3 if
    % dx=0.3125 cm results
    HRR3 = importdata([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_p3125cm_hrr.csv'],',',2);
    DEV3 = importdata([outdir,fuel_name{i_fuel},'_XO2_ramp_dx_p3125cm_devc.csv'],',',2);

    Time_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'Time')));
    XO2_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'"XO2"')));
    q_R_FDS_3 = 0.5*( DEV3.data(:,find(strcmp(DEV3.colheaders,'"qrad1"'))) + DEV3.data(:,find(strcmp(DEV3.colheaders,'"qrad2"'))) );
    HRR_FDS_3 = HRR3.data(:,find(strcmp(HRR3.colheaders,'HRR')));
    MLR_FDS_3 = HRR3.data(:,find(strcmp(HRR3.colheaders,'MLR_FUEL')));
    eta_FDS_3 = HRR_FDS_3./(MLR_FDS_3*fuel_hoc(i_fuel));
    CHI_R_FDS_3 = DEV3.data(:,find(strcmp(DEV3.colheaders,'"CHI_R"')));
    QR_FDS_3 = abs(HRR3.data(:,find(strcmp(HRR3.colheaders,'Q_RADI'))));
    GLOB_CHI_R_3 = QR_FDS_3./HRR_FDS_3;

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
    Lf_tmp = Lf_FDS_3;
    for n=1:length(Time_FDS_3)
        indx_range = [find(Time_FDS_3>(Time_FDS_3(n)-Lf_dt),1):n];
        Lf_FDS_3(n) = mean(Lf_tmp(indx_range));
    end
end % case 3 if

if is_case_4_run % case 4 if
    % 2 step dx=1.25 cm results
    HRR4 = importdata([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_1p25cm_hrr.csv'],',',2);
    DEV4 = importdata([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_1p25cm_devc.csv'],',',2);

    Time_FDS_4 = DEV4.data(:,find(strcmp(DEV4.colheaders,'Time')));
    XO2_FDS_4 = DEV4.data(:,find(strcmp(DEV4.colheaders,'"XO2"')));
    q_R_FDS_4 = 0.5*( DEV4.data(:,find(strcmp(DEV4.colheaders,'"qrad1"'))) + DEV4.data(:,find(strcmp(DEV4.colheaders,'"qrad2"'))) );
    HRR_FDS_4 = HRR4.data(:,find(strcmp(HRR4.colheaders,'HRR')));
    MLR_FDS_4 = HRR4.data(:,find(strcmp(HRR4.colheaders,'MLR_FUEL')));
    eta_FDS_4 = HRR_FDS_4./(MLR_FDS_4*fuel_hoc(i_fuel));
    CHI_R_FDS_4 = DEV4.data(:,find(strcmp(DEV4.colheaders,'"CHI_R"')));
    QR_FDS_4 = abs(HRR4.data(:,find(strcmp(HRR4.colheaders,'Q_RADI'))));
    GLOB_CHI_R_4 = QR_FDS_4./HRR_FDS_4;

    % flame height
    colLf01 = find(strcmp(DEV4.colheaders,'"Lf-01"'));
    colLf50 = find(strcmp(DEV4.colheaders,'"Lf-50"'));
    Lf_FDS_4 = zeros(1,length(Time_FDS_4));
    for n=1:length(Time_FDS_4)
        hrrpul = DEV4.data(n,colLf01:colLf50);
        q_total = sum(hrrpul);
        q_sum = 0.;
        for i=1:Lf_pts
            q_sum = q_sum+hrrpul(i);
            if q_sum>Lf_fac*q_total
                Lf_FDS_4(n) = Lf_z(i);
                break;
            end
        end
    end
    Lf_tmp = Lf_FDS_4;
    for n=1:length(Time_FDS_4)
        indx_range = [find(Time_FDS_4>(Time_FDS_4(n)-Lf_dt),1):n];
        Lf_FDS_4(n) = mean(Lf_tmp(indx_range));
    end
end % case 4 if

if is_case_5_run % case 5 if
    % 2 step dx=0.625 cm results
    HRR5 = importdata([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_p625cm_hrr.csv'],',',2);
    DEV5 = importdata([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_p625cm_devc.csv'],',',2);

    Time_FDS_5 = DEV5.data(:,find(strcmp(DEV5.colheaders,'Time')));
    XO2_FDS_5 = DEV5.data(:,find(strcmp(DEV5.colheaders,'"XO2"')));
    q_R_FDS_5 = 0.5*( DEV5.data(:,find(strcmp(DEV5.colheaders,'"qrad1"'))) + DEV5.data(:,find(strcmp(DEV5.colheaders,'"qrad2"'))) );
    HRR_FDS_5 = HRR5.data(:,find(strcmp(HRR5.colheaders,'HRR')));
    MLR_FDS_5 = HRR5.data(:,find(strcmp(HRR5.colheaders,'MLR_FUEL')));
    eta_FDS_5 = HRR_FDS_5./(MLR_FDS_5*fuel_hoc(i_fuel));
    CHI_R_FDS_5 = DEV5.data(:,find(strcmp(DEV5.colheaders,'"CHI_R"')));
    QR_FDS_5 = abs(HRR5.data(:,find(strcmp(HRR5.colheaders,'Q_RADI'))));
    GLOB_CHI_R_5 = QR_FDS_5./HRR_FDS_5;

    % flame height
    colLf01 = find(strcmp(DEV5.colheaders,'"Lf-01"'));
    colLf50 = find(strcmp(DEV5.colheaders,'"Lf-50"'));
    Lf_FDS_5 = zeros(1,length(Time_FDS_5));
    for n=1:length(Time_FDS_5)
        hrrpul = DEV5.data(n,colLf01:colLf50);
        q_total = sum(hrrpul);
        q_sum = 0.;
        for i=1:Lf_pts
            q_sum = q_sum+hrrpul(i);
            if q_sum>Lf_fac*q_total
                Lf_FDS_5(n) = Lf_z(i);
                break;
            end
        end
    end
    Lf_tmp = Lf_FDS_5;
    for n=1:length(Time_FDS_5)
        indx_range = [find(Time_FDS_5>(Time_FDS_5(n)-Lf_dt),1):n];
        Lf_FDS_5(n) = mean(Lf_tmp(indx_range));
    end
end % case 5 if

if is_case_6_run % case 6 if
    % 2 step dx=0.3125 cm results
    HRR6 = importdata([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_p3125cm_hrr.csv'],',',2);
    DEV6 = importdata([outdir,fuel_name{i_fuel},'_2step_XO2_ramp_dx_p3125cm_devc.csv'],',',2);

    Time_FDS_6 = DEV6.data(:,find(strcmp(DEV6.colheaders,'Time')));
    XO2_FDS_6 = DEV6.data(:,find(strcmp(DEV6.colheaders,'"XO2"')));
    q_R_FDS_6 = 0.5*( DEV6.data(:,find(strcmp(DEV6.colheaders,'"qrad1"'))) + DEV6.data(:,find(strcmp(DEV6.colheaders,'"qrad2"'))) );
    HRR_FDS_6 = HRR6.data(:,find(strcmp(HRR6.colheaders,'HRR')));
    MLR_FDS_6 = HRR6.data(:,find(strcmp(HRR6.colheaders,'MLR_FUEL')));
    eta_FDS_6 = HRR_FDS_6./(MLR_FDS_6*fuel_hoc(i_fuel));
    CHI_R_FDS_6 = DEV6.data(:,find(strcmp(DEV6.colheaders,'"CHI_R"')));
    QR_FDS_6 = abs(HRR6.data(:,find(strcmp(HRR6.colheaders,'Q_RADI'))));
    GLOB_CHI_R_6 = QR_FDS_6./HRR_FDS_6;

    % flame height
    colLf01 = find(strcmp(DEV6.colheaders,'"Lf-01"'));
    colLf50 = find(strcmp(DEV6.colheaders,'"Lf-50"'));
    Lf_FDS_6 = zeros(1,length(Time_FDS_6));
    for n=1:length(Time_FDS_6)
        hrrpul = DEV6.data(n,colLf01:colLf50);
        q_total = sum(hrrpul);
        q_sum = 0.;
        for i=1:Lf_pts
            q_sum = q_sum+hrrpul(i);
            if q_sum>Lf_fac*q_total
                Lf_FDS_6(n) = Lf_z(i);
                break;
            end
        end
    end
    Lf_tmp = Lf_FDS_6;
    for n=1:length(Time_FDS_6)
        indx_range = [find(Time_FDS_6>(Time_FDS_6(n)-Lf_dt),1):n];
        Lf_FDS_6(n) = mean(Lf_tmp(indx_range));
    end
end % case 6 if

% compute CHI_R ramp
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

X_O2_max = 2.077972E-01;
X_O2_min = 0.1;
ramp_time = @(x) (x-X_O2_max)*60./(X_O2_min-X_O2_max); % this comes from linear XO2 ramp (see input file)
H(1)=plot(XO2,Chi_R); hold on
switch i_fuel
    case 1 % methane Chi_r ramp
        pw_ramp_time = [0,37,45,60];
        pw_ramp_XO2  = [0.21,0.142,0.123,0];
        pw_ramp_chir = [.245,.125,0,0];
    case 2
        pw_ramp_time = [0,37,45,60];
        pw_ramp_XO2  = [0.21,.17,0.133,0.129,0];
        pw_ramp_chir = [.33,.245,.125,0,0];
end
H(2)=plot(pw_ramp_XO2,pw_ramp_chir);

if is_case_1_run | is_case_2_run | is_case_3_run | is_case_4_run | is_case_5_run | is_case_6_run
    is_case_run = 1;
else
    is_case_run = 0;
end

if is_case_run
    if is_case_1_run; H(3)=plot(XO2_FDS_1,CHI_R_FDS_1,line_fmt{1}); end
    if is_case_2_run; H(4)=plot(XO2_FDS_2,CHI_R_FDS_2,line_fmt{2}); end
    if is_case_3_run; H(5)=plot(XO2_FDS_3,CHI_R_FDS_3,line_fmt{3}); end
    if is_case_4_run; H(6)=plot(XO2_FDS_4,CHI_R_FDS_4,line_fmt{4}); end
    if is_case_5_run; H(7)=plot(XO2_FDS_5,CHI_R_FDS_5,line_fmt{5}); end
    if is_case_6_run; H(8)=plot(XO2_FDS_6,CHI_R_FDS_6,line_fmt{6}); end
    axis([0.12 0.22 0 0.35 ])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('\chi_R','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    text(5,0.27,'Radiative Fraction Ramp','FontName',Font_Name,'FontSize',Title_Font_Size)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    lh=legend(H,'Exp','Linear Ramp',key_fmt{:},'Location','SouthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_Chi_r_ramp_check']);
end

clear H
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

steel_blue = [0 0.447 0.741];
H(1)=plot(XO2,eta,'.','MarkerSize',10); hold on
set(H(1),'Color',steel_blue)
subr= 1:10:length(XO2);
h=errorbar(XO2(subr),eta(subr),-S_eta(subr),S_eta(subr),-S_XO2(subr),S_XO2(subr),'.','MarkerSize',10); hold on
set(h,'Color',steel_blue)
if is_case_run
    if is_case_1_run; H(2)=plot(XO2_FDS_1,eta_FDS_1,line_fmt{1}); end
    if is_case_2_run; H(3)=plot(XO2_FDS_2,eta_FDS_2,line_fmt{2}); end
    if is_case_3_run; H(4)=plot(XO2_FDS_3,eta_FDS_3,line_fmt{3}); end
    if is_case_4_run; H(5)=plot(XO2_FDS_4,eta_FDS_4,line_fmt{4}); end
    if is_case_5_run; H(6)=plot(XO2_FDS_5,eta_FDS_5,line_fmt{5}); end
    if is_case_6_run; H(7)=plot(XO2_FDS_6,eta_FDS_6,line_fmt{6}); end
    axis([0.09 0.21 0 1.2 ])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('\eta','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    text(0.095,1.12,[Fuel_name{i_fuel},' Combustion Efficiency'],'FontName',Font_Name,'FontSize',Title_Font_Size)

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

    lh=legend(H,'Exp',key_fmt{:},'Location','SouthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_eta']);
end

% double check N2 ramp

% intended ramp:
Time_ramp = [0,60];
XO2_ramp = [X_O2_max,X_O2_min];
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(Time_ramp,XO2_ramp,'--o','MarkerSize',10); hold on
if is_case_run
    if is_case_1_run; H(2)=plot(Time_FDS_1,XO2_FDS_1,line_fmt{1}); end
    if is_case_2_run; H(3)=plot(Time_FDS_2,XO2_FDS_2,line_fmt{2}); end
    if is_case_3_run; H(4)=plot(Time_FDS_3,XO2_FDS_3,line_fmt{3}); end
    if is_case_4_run; H(5)=plot(Time_FDS_4,XO2_FDS_4,line_fmt{4}); end
    if is_case_5_run; H(6)=plot(Time_FDS_5,XO2_FDS_5,line_fmt{5}); end
    if is_case_6_run; H(7)=plot(Time_FDS_6,XO2_FDS_6,line_fmt{6}); end
    axis([0 60 0.1 0.21])
    xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    lh=legend(H,'Exp',key_fmt{:},'Location','SouthWest');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_N2_ramp_check']);
end

% plot heat flux
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(XO2,q_R,'o'); hold on
if is_case_run
    if is_case_1_run; H(2)=plot(XO2_FDS_1,q_R_FDS_1,line_fmt{1}); end
    if is_case_2_run; H(3)=plot(XO2_FDS_2,q_R_FDS_2,line_fmt{2}); end
    if is_case_3_run; H(4)=plot(XO2_FDS_3,q_R_FDS_3,line_fmt{3}); end
    if is_case_4_run; H(5)=plot(XO2_FDS_4,q_R_FDS_4,line_fmt{4}); end
    if is_case_5_run; H(6)=plot(XO2_FDS_5,q_R_FDS_5,line_fmt{5}); end
    if is_case_6_run; H(7)=plot(XO2_FDS_6,q_R_FDS_6,line_fmt{6}); end
    axis([0.12 0.21 0 1.5])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('{\itq"}_R (kW/m^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    lh=legend(H,'Exp',key_fmt{:},'Location','SouthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_rad_heat_flux']);
end

% plot flame height
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(XO2_Lf,Lf,'o'); hold on
subr=1:2:length(Lf);
h=errorbar(XO2_Lf(subr),Lf(subr),-S_Lf(subr),S_Lf(subr),'o'); hold on
set(h,'Color',steel_blue)
if is_case_run
    if is_case_1_run; H(2)=plot(XO2_FDS_1,Lf_FDS_1,line_fmt{1}); end
    if is_case_2_run; H(3)=plot(XO2_FDS_2,Lf_FDS_2,line_fmt{2}); end
    if is_case_3_run; H(4)=plot(XO2_FDS_3,Lf_FDS_3,line_fmt{3}); end
    if is_case_4_run; H(5)=plot(XO2_FDS_4,Lf_FDS_4,line_fmt{4}); end
    if is_case_5_run; H(6)=plot(XO2_FDS_5,Lf_FDS_5,line_fmt{5}); end
    if is_case_6_run; H(7)=plot(XO2_FDS_6,Lf_FDS_6,line_fmt{6}); end
    axis([0.14 0.22 0.4 1.00])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('{\itL}_f (m)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    lh=legend(H,'Exp',key_fmt{:},'Location','NorthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_flame_height']);
end

% plot actual global radiative fraction
figure
clear H
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
H(1)=plot(XO2,Chi_R); hold on
H(2)=plot(pw_ramp_XO2,pw_ramp_chir);
if is_case_run
    if is_case_1_run; H(3)=plot(XO2_FDS_1,GLOB_CHI_R_1,line_fmt{1}); end
    if is_case_2_run; H(4)=plot(XO2_FDS_2,GLOB_CHI_R_2,line_fmt{2}); end
    if is_case_3_run; H(5)=plot(XO2_FDS_3,GLOB_CHI_R_3,line_fmt{3}); end
    if is_case_4_run; H(6)=plot(XO2_FDS_4,GLOB_CHI_R_4,line_fmt{4}); end
    if is_case_5_run; H(7)=plot(XO2_FDS_5,GLOB_CHI_R_5,line_fmt{5}); end
    if is_case_6_run; H(8)=plot(XO2_FDS_6,GLOB_CHI_R_6,line_fmt{6}); end
    axis([0.12 0.22 0 0.35 ])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('\chi_R','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    text(5,0.27,'Radiative Fraction Ramp','FontName',Font_Name,'FontSize',Title_Font_Size)
    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    lh=legend(H,'Exp','Linear Ramp',key_fmt{:},'Location','SouthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_global_Chi_R']);
end

