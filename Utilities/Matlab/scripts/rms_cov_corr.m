% Floyd
% 5-23-2014
% rms_cov_corr.m

close all
clear all

addpath('../../Verification/Controls');

% load experimental data and FDS prediction
M = importdata('rms_cov_corr_devc.csv',',',2);
fds_data = M.data;

startrow=100/0.02+1;
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

xcalc(1)=60;
xcalc(2)=1000;
ycalc(1)=urms;
ycalc(2)=urms;
maxval=max(fds_data(:,4));
maxval=ceil(maxval*10)/10;

h=plot(xcalc,ycalc,'r-',fds_data(:,1),fds_data(:,4),'k-','LineWidth',1.5);
hold on
axis([0 1000 0 maxval])
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
xlabel('Time(s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it u} rms (m/s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
legend('Analytic','FDS','Location','SouthWest')

svn_file = '../../Verification/Controls/rms_cov_corr_svn.txt';

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
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
plotname = ['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/rms_cov_corr_rms'];
print(gcf,'-dpdf',plotname);

hold off
clear h

ycalc(1)=uwcov;
ycalc(2)=uwcov;
maxval=max(fds_data(:,5));
maxval=ceil(maxval*10)/10;

h=plot(xcalc,ycalc,'r-',fds_data(:,1),fds_data(:,5),'k-','LineWidth',1.5);
hold on
axis([0 1000 0 maxval])
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
xlabel('Time(s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it uw} covariance (m^2/s^2)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
legend('Analytic','FDS','Location','NorthEast')

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
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
plotname = ['../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/rms_cov_corr_cov'];
print(gcf,'-dpdf',plotname);

hold off
clear h

ycalc(1)=uwcorr;
ycalc(2)=uwcorr;
maxval= max(fds_data(:,6));
maxval=ceil(maxval*10)/10;

h=plot(xcalc,ycalc,'r-',fds_data(:,1),fds_data(:,6),'k-','LineWidth',1.5);
hold on
axis([0 1000 0 maxval])
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
xlabel('Time(s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('{\it uw} cross correlation','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
legend('Analytic','FDS','Location','NorthEast')

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
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
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
