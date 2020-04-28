% Floyd
% 5-23-2014
% rms_cov_corr.m

close all
clear all

plot_style

datadir='../../Verification/Controls/';

% load experimental data and FDS prediction

filename = [datadir,'rms_cov_corr_devc.csv'];

if ~exist(filename) % skip_case_if

    display(['Error: File ' filename ' does not exist. Skipping case.'])

else

    fds_data = csvread([datadir,'rms_cov_corr_devc.csv'],2);

    startrow=500/0.02+1;
    endrow=size(fds_data,1);

    umean=mean(fds_data(startrow:endrow,2));
    wmean=mean(fds_data(startrow:endrow,3));

    udiff=fds_data(startrow:endrow,2)-umean;
    wdiff=fds_data(startrow:endrow,3)-wmean;

    udiff2=udiff;
    wdiff2=wdiff;
    uwcova=udiff;
    for i=1:endrow-startrow+1
        udiff2(i)=udiff(i)*udiff(i);
        wdiff2(i)=wdiff(i)*wdiff(i);
        uwcova(i)=udiff(i)*wdiff(i);
        drawnow;
    end

    urms=sqrt(mean(udiff2));
    wrms=sqrt(mean(wdiff2));
    uwcov=mean(uwcova);
    uwcorr=uwcov/urms/wrms;

    urms_fds=fds_data(endrow,4);
    uwcov_fds=fds_data(endrow,5);
    uwcorr_fds=fds_data(endrow,6);

    xcalc(1)=0;
    xcalc(2)=1000;
    ycalc(1)=urms;
    ycalc(2)=urms;

    maxval=ceil(2*urms*100)/100;

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    h=plot(xcalc,ycalc,'r-',fds_data(:,1),fds_data(:,4),'k-');

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    xlabel('Time(s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
    ylabel('{\it u} rms (m/s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
    legend('Analytic','FDS','Location','East')

    % add Git revision if file is available
    git_file = '../../Verification/Controls/rms_cov_corr_git.txt';
    addverstr(gca,git_file,'linear')

    % print to pdf
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    plotname = ['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/rms_cov_corr_rms'];
    print(gcf,'-dpdf',plotname);

    clear h

    ycalc(1)=uwcov;
    ycalc(2)=uwcov;

    maxval=ceil(2*uwcov*100)/100;

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    h=plot(xcalc,ycalc,'r-',fds_data(:,1),fds_data(:,5),'k-');

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    xlabel('Time(s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
    ylabel('{\it uw} covariance (m^2/s^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
    legend('Analytic','FDS','Location','East')

    % add Git revision if file is available
    git_file = '../../Verification/Controls/rms_cov_corr_git.txt';
    addverstr(gca,git_file,'linear')

    % print to pdf
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    plotname = ['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/rms_cov_corr_cov'];
    print(gcf,'-dpdf',plotname);

    clear h

    ycalc(1)=uwcorr;
    ycalc(2)=uwcorr;

    maxval=ceil(2*uwcorr*100)/100;

    figure
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    h=plot(xcalc,ycalc,'r-',fds_data(:,1),fds_data(:,6),'k-');

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    xlabel('Time(s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
    ylabel('{\it uw} cross correlation','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'FontName',Font_Name)
    legend('Analytic','FDS','Location','East')

    % add Git revision if file is available
    git_file = '../../Verification/Controls/rms_cov_corr_git.txt';
    addverstr(gca,git_file,'linear')

    % print to pdf
    set(gcf,'Visible',Figure_Visibility);
    set(gcf,'Units',Paper_Units);
    set(gcf,'PaperUnits',Paper_Units);
    set(gcf,'PaperSize',[Paper_Width Paper_Height]);
    set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
    plotname = ['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/rms_cov_corr_corr'];
    print(gcf,'-dpdf',plotname);

    clear h

    % check errors
    if abs((urms-urms_fds)/urms) > 0.001
       display(['Matlab Warning: urms in rms_cov_corr is out of tolerance.'])
    end
    if abs((uwcov-uwcov_fds)/urms) > 0.001
       display(['Matlab Warning: uwcov in rms_cov_corr is out of tolerance.'])
    end
    if abs((uwcorr-uwcorr_fds)/urms) > 0.001
       display(['Matlab Warning: uwcorr in rms_cov_corr is out of tolerance.'])
    end

end % skip_case_if
