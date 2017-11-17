% McDermott
% 6 Sep 2017
% umd_line_burner_3.m
%
% Plot global radiative fraction and combustion efficiency for UMD Line Burner cases.

close all
clear all

plot_style

expdir = '../../../exp/Submodules/macfp-db/Extinction/UMD_Line_Burner/Experimental_Data/';
outdir = '../../../out/UMD_Line_Burner/FDS_Output_Files/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/UMD_Line_Burner/';

exp_fname    = {'CH4_A_Data.csv','C3H8_A_Data.csv'};
fuel_name    = {'methane','propane'};
Fuel_name    = {'Methane','Propane'};
fuel_hoc     = [50010.3475,46334.6246]; % from .out file
git_tag_ext  = '_XO2_ramp_dx_1p25cm_git.txt';

% line_fmt = {'ro-','mo-','bo-','r^-','m^-','b^-'};
% key_fmt  = {'1step, dx=1.25cm','1step, dx=0.625cm','1step, dx=0.3125cm','2step, dx=1.25cm','2step, dx=0.625cm','2step, dx=0.3125cm'};

% case_name = {'_XO2_ramp_dx_1p25cm',...
%              '_XO2_ramp_dx_p625cm',...
%              '_XO2_ramp_dx_p3125cm',...
%              '_2step_XO2_ramp_dx_1p25cm',...
%              '_2step_XO2_ramp_dx_p625cm',...
%              '_2step_XO2_ramp_dx_p3125cm'};

% omit high res cases ---------------------------
line_fmt = {'ro-','mo-','r^-','m^-'};
key_fmt  = {'1step, dx=1.25cm','1step, dx=0.625cm','2step, dx=1.25cm','2step, dx=0.625cm',};

case_name = {'_XO2_ramp_dx_1p25cm',...
             '_XO2_ramp_dx_p625cm',...
             '_2step_XO2_ramp_dx_1p25cm',...
             '_2step_XO2_ramp_dx_p625cm'};
% -----------------------------------------------

for i_fuel = 1:2 % i_fuel_loop

    % experimental results
    EXP = importdata([expdir,exp_fname{i_fuel}],',',1);
    XO2   = EXP.data(:,find(strcmp(EXP.colheaders,'XO2')));
    q_R   = EXP.data(:,find(strcmp(EXP.colheaders,'q_R')));
    Chi_R = EXP.data(:,find(strcmp(EXP.colheaders,'Chi_R')));
    S_XO2 = EXP.data(:,find(strcmp(EXP.colheaders,'S_XO2'))); % uncertainty
    eta   = EXP.data(:,find(strcmp(EXP.colheaders,'eta')));
    S_eta = EXP.data(:,find(strcmp(EXP.colheaders,'S_eta')));

    figure(1)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
    Hfig1(1)=plot(XO2,Chi_R,'k.','MarkerSize',15); hold on

    figure(2)
    Hfig2(1)=plot(XO2,eta,'k.','MarkerSize',15); hold on
    subr= 1:10:length(XO2);
    herr=errorbar(XO2(subr),eta(subr),-S_eta(subr),S_eta(subr),-S_XO2(subr),S_XO2(subr),'.','MarkerSize',10);
    set(herr,'Color','k')

    for j_case=1:length(case_name)

        if ~exist([outdir,fuel_name{i_fuel},case_name{j_case},'_hrr.csv'])
            disp('Matlab Warning: ',[outdir,fuel_name{i_fuel},case_name{j_case},'_hrr.csv'],' does not exist.')
            continue
        end
        if ~exist([outdir,fuel_name{i_fuel},case_name{j_case},'_devc.csv'])
            disp('Matlab Warning: ',[outdir,fuel_name{i_fuel},case_name{j_case},'_devc.csv'],' does not exist.')
            continue
        end

        HRR = importdata([outdir,fuel_name{i_fuel},case_name{j_case},'_hrr.csv'],',',2);
        DEV = importdata([outdir,fuel_name{i_fuel},case_name{j_case},'_devc.csv'],',',2);

        %Time_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'Time')));
        XO2_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'"XO2"')));
        %q_R_FDS = 0.5*( DEV.data(:,find(strcmp(DEV.colheaders,'"qrad1"'))) + DEV.data(:,find(strcmp(DEV.colheaders,'"qrad2"'))) );
        HRR_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'HRR')));
        MLR_FDS = HRR.data(:,find(strcmp(HRR.colheaders,'MLR_FUEL')));
        eta_FDS = HRR_FDS./(MLR_FDS*fuel_hoc(i_fuel));
        %CHI_R_FDS = DEV.data(:,find(strcmp(DEV.colheaders,'"CHI_R"')));
        QR_FDS = abs(HRR.data(:,find(strcmp(HRR.colheaders,'Q_RADI'))));

        % clip low hrr to zero to avoid spuriously high results (converts to Inf)
        for k=1:length(HRR_FDS)
            if HRR_FDS(k)<5
                HRR_FDS(k)=0;
            end
        end
        GLOB_CHI_R = QR_FDS./HRR_FDS;

        figure(1)
        Hfig1(1+j_case)=plot(XO2_FDS,GLOB_CHI_R,line_fmt{j_case});

        figure(2)
        Hfig2(1+j_case)=plot(XO2_FDS,eta_FDS,line_fmt{j_case});

    end

    figure(1)
    hold off
    axis([0.12 0.22 0 0.35 ])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('\chi_R','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    text(0.125,0.325,[Fuel_name{i_fuel},' Global Radiative Fraction'],'FontName',Font_Name,'FontSize',Title_Font_Size)

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    lh=legend(Hfig1,'Exp',key_fmt{:},'Location','SouthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)
    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_global_Chi_R']);

    figure(2)
    hold off
    axis([0.09 0.21 0 1.2 ])
    xlabel('O2 (vol frac)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    ylabel('\eta','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
    text(0.095,1.12,[Fuel_name{i_fuel},' Combustion Efficiency'],'FontName',Font_Name,'FontSize',Title_Font_Size)

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)

    lh=legend(Hfig2,'Exp',key_fmt{:},'Location','SouthEast');
    set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

    git_file=[outdir,fuel_name{i_fuel},git_tag_ext];
    addverstr(gca,git_file,'linear');

    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    print(gcf,'-dpdf',[pltdir,fuel_name{i_fuel},'_eta']);

end % i_fuel_loop

%%%%%%%%% SAVE %%%%%%%%%
% % compute CHI_R ramp
% figure
% set(gca,'Units',Plot_Units)
% set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% X_O2_max = 2.077972E-01;
% X_O2_min = 0.1;
% ramp_time = @(x) (x-X_O2_max)*60./(X_O2_min-X_O2_max); % this comes from linear XO2 ramp (see input file)
% H(1)=plot(XO2,Chi_R); hold on
% switch i_fuel
%     case 1 % methane Chi_r ramp
%         pw_ramp_time = [0,37,45,60];
%         pw_ramp_XO2  = [0.21,0.142,0.123,0];
%         pw_ramp_chir = [.245,.125,0,0];
%     case 2
%         pw_ramp_time = [0,37,45,60];
%         pw_ramp_XO2  = [0.21,.17,0.133,0.129,0];
%         pw_ramp_chir = [.33,.245,.125,0,0];
% end
% H(2)=plot(pw_ramp_XO2,pw_ramp_chir);


