% McDermott
% 5-14-2009
% poiseuille_convergence.m

close all
clear all

dpdx = -1;
L = 1;
N = [8,16,32,64];

dir = '../../Validation/Moody_Chart/FDS_Output_Files/';

[f(1),Re(1)] = friction_factor_calc(dpdx,L,[dir,'poiseuille_N8_mu025_devc.csv']);
[f(2),Re(2)] = friction_factor_calc(dpdx,L,[dir,'poiseuille_N16_mu025_devc.csv']);
[f(3),Re(3)] = friction_factor_calc(dpdx,L,[dir,'poiseuille_N32_mu025_devc.csv']);
[f(4),Re(4)] = friction_factor_calc(dpdx,L,[dir,'poiseuille_N64_mu025_devc.csv']);

% plot convergence for Poiseuille flow (mu = 0.025)

dz = L./N;
error = abs(f-24./Re);
H(1)=loglog(dz,error,'b*-','Linewidth',1.); hold on
H(2)=loglog(dz,.05*dz,'k--','Linewidth',1.);
H(3)=loglog(dz,.4*dz.^2,'k-','Linewidth',1.);

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0.01 0.2 0.00005 0.01])

xlabel('Grid Spacing, {\it \deltaz} (m)','Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('Friction Factor Error')
h = legend(H,'FDS','O({\it \deltaz})','O({\it \deltaz}^2)','Location','Southeast');
set(h,'Interpreter',Font_Interpreter)

% add SVN if file is available

SVN_Filename = [dir,'poiseuille_N8_mu025_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = 10^( log10(x_lim(1))+ SVN_Scale_X*( log10(x_lim(2)) - log10(x_lim(1)) ) );
    Y_SVN_Position = 10^( log10(y_lim(1))+ SVN_Scale_Y*( log10(y_lim(2)) - log10(y_lim(1)) ) );
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/poiseuille_convergence')
