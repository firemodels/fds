% Abbas
% 7-3-12
% yplus.m

close all
clear all

dir = '../../Verification/Turbulence/';

skip_case = 0;
if ~exist([dir,'yplus_8_devc.csv'])
    display(['Error: File ' [dir,'yplus_8_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([dir,'yplus_16_devc.csv'])
    display(['Error: File ' [dir,'yplus_16_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([dir,'yplus_32_devc.csv'])
    display(['Error: File ' [dir,'yplus_32_devc.csv'] ' does not exist. Skipping case.'])
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

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(yp_exact,yp_exact,'ro'); hold on
plot(yp_exact,yp,'b--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it y}^+ specified','Interpreter',Font_Interpreter,'Fontname',Font_Name),
ylabel('{\it y}^+ predicted','Interpreter',Font_Interpreter,'Fontname',Font_Name),
leg = legend('exact','FDS');
set(leg,'Location','SouthEast')
set(leg,'Interpreter',Font_Interpreter)
set(leg,'FontName',Font_Name)
set(leg,'FontSize',Key_Font_Size)

% add git version if file is available

git_file = [dir,'yplus_8_git.txt'];
addverstr(gca,git_file,'linear')

% print pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/yplus')

