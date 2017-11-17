% McDermott
% 5-14-2009
% poiseuille_convergence.m

close all
clear all

plot_style

dpdx = -1;
L = 1;
N = [8,16,32,64];

outdir = '../../../out/Moody_Chart/FDS_Output_Files/';

[f(1),Re(1)] = friction_factor_calc(dpdx,L,[outdir,'poiseuille_N8_mu025_devc.csv']);
[f(2),Re(2)] = friction_factor_calc(dpdx,L,[outdir,'poiseuille_N16_mu025_devc.csv']);
[f(3),Re(3)] = friction_factor_calc(dpdx,L,[outdir,'poiseuille_N32_mu025_devc.csv']);
[f(4),Re(4)] = friction_factor_calc(dpdx,L,[outdir,'poiseuille_N64_mu025_devc.csv']);

% plot convergence for Poiseuille flow (mu = 0.025)

dz = L./N;
error = abs(f-24./Re);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=loglog(dz,error,'b*-','Linewidth',1.); hold on
H(2)=loglog(dz,.05*dz,'k--','Linewidth',1.);
H(3)=loglog(dz,.4*dz.^2,'k-','Linewidth',1.);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0.01 0.2 0.00005 0.01])

xlabel('Grid Spacing, {\it\deltaz} (m)','Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('Friction Factor Error')
h = legend(H,'FDS','{\itO}({\it\deltaz})','{\itO}({\it\deltaz}^2)','Location','Southeast');
set(h,'Interpreter',Font_Interpreter)

% add Git revision if file is available

Git_Filename = [outdir,'poiseuille_N8_mu025_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/poiseuille_convergence')
