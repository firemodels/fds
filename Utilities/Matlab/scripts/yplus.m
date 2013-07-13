% Abbas
% 7-3-12
% yplus.m

close all
clear all

dir = '../../Verification/Turbulence/';

skip_case = 0;
if ~exist([dir,'yplus_8_devc.csv'])
    display(['Error: Files ' [dir,'yplus_8_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([dir,'yplus_16_devc.csv'])
    display(['Error: Files ' [dir,'yplus_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([dir,'yplus_32_devc.csv'])
    display(['Error: Files ' [dir,'yplus_32_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if skip_case
    return
end

M_8  = importdata([dir,'yplus_8_devc.csv'],',',2);
M_16 = importdata([dir,'yplus_16_devc.csv'],',',2);
M_32 = importdata([dir,'yplus_32_devc.csv'],',',2);

yp = [M_8.data(end,2),M_16.data(end,2),M_32.data(end,2)];

n = [8 16 32];
y = 0.5*(1./n);
mu = 1;
rho = 1.199;
tau_w = 0.5;
u_tau = sqrt(tau_w/rho);
d_nu = (mu/rho)/u_tau;
yp_exact = y/d_nu;

plot(yp_exact,yp_exact,'ro'); hold on
plot(yp_exact,yp,'b--')

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)

xlabel('{\it y}^+ specified','Interpreter',Font_Interpreter,'Fontname',Font_Name),
ylabel('{\it y}^+ predicted','Interpreter',Font_Interpreter,'Fontname',Font_Name),
leg = legend('exact','FDS');
set(leg,'Location','SouthEast')
set(leg,'Interpreter',Font_Interpreter)
set(leg,'Fontname',Font_Name)

% add SVN if file is available

SVN_Filename = [dir,'yplus_8_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/yplus')

