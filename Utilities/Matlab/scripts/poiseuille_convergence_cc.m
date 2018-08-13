% CC script that follows McDermott poiseuille_convergence.m

close all
clear all

plot_style

dpdx = -1;
L    = 1;
N    = [10,20,40,80];

outdir = '../../Verification/Complex_Geometry/';

% Aligned case theta=0:
[f(1),Re(1)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N10a_theta0_devc.csv']);
[f(2),Re(2)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N20a_theta0_devc.csv']);
[f(3),Re(3)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N40a_theta0_devc.csv']);
[f(4),Re(4)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N80a_theta0_devc.csv']);

[f2(1),Re2(1)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N10nah_theta0_devc.csv']);
[f2(2),Re2(2)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N20nah_theta0_devc.csv']);
[f2(3),Re2(3)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N40nah_theta0_devc.csv']);
[f2(4),Re2(4)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N80nah_theta0_devc.csv']);

% plot convergence for Poiseuille flow aligned case theta=0 (mu = 0.025)

dz = L./N;
error = abs(f-24./Re);
error2= abs(f2-24./Re2);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=loglog(dz,error,'b*-','Linewidth',1.); hold on
H(2)=loglog(dz,error2,'rx-','Linewidth',1.); hold on
H(3)=loglog(dz,.12*dz,'k--','Linewidth',1.);
H(4)=loglog(dz,.4*dz.^2,'k-','Linewidth',1.);

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0.01 0.2 0.0000005 0.01])

xlabel('Grid Spacing, {\it\deltaz} (m)','Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('Friction Factor Error')
h = legend(H,'FDS, h=0','FDS, h=\Deltaz/3','{\itO}({\it\deltaz})','{\itO}({\it\deltaz}^2)','Location','Southeast');
set(h,'Interpreter',Font_Interpreter)

% add Git revision if file is available

Git_Filename = [outdir,'geom_poiseuille_N10a_theta0_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/geom_poiseuille_convergence_theta0a')


% NOT aligned case theta=0

clear f Re f2 Re2 H

[f(1),Re(1)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N10na_theta0_devc.csv']);
[f(2),Re(2)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N20na_theta0_devc.csv']);
[f(3),Re(3)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N40na_theta0_devc.csv']);
[f(4),Re(4)] = friction_factor_calc(dpdx,L,[outdir,'geom_poiseuille_N80na_theta0_devc.csv']);

% plot convergence for Poiseuille flow not aligned case theta=0 (mu = 0.025)

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
axis([0.01 0.2 0.0000005 0.01])

xlabel('Grid Spacing, {\it\deltaz} (m)','Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('Friction Factor Error')
h = legend(H,'FDS, h=\Deltaz_{10}/11','{\itO}({\it\deltaz})','{\itO}({\it\deltaz}^2)','Location','Southeast');
set(h,'Interpreter',Font_Interpreter)

% add Git revision if file is available

Git_Filename = [outdir,'geom_poiseuille_N10na_theta0_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/geom_poiseuille_convergence_theta0na')


