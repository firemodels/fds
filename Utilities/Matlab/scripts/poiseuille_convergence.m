% McDermott
% 5-14-2009
% poiseuille_convergence.m

close all
clear all

dpdx = -1;
L = 1;
N = [8,16,32,64];

[f(1),Re(1)] = friction_factor_calc(dpdx,L,'../../../Verification/Flowfields/poiseuille_N8_mu025_devc.csv');
[f(2),Re(2)] = friction_factor_calc(dpdx,L,'../../../Verification/Flowfields/poiseuille_N16_mu025_devc.csv');
[f(3),Re(3)] = friction_factor_calc(dpdx,L,'../../../Verification/Flowfields/poiseuille_N32_mu025_devc.csv');
[f(4),Re(4)] = friction_factor_calc(dpdx,L,'../../../Verification/Flowfields/poiseuille_N64_mu025_devc.csv');

% plot convergence for Poiseuille flow (mu = 0.025)

dz = L./N;
error = abs(f-24./Re);
H(1)=loglog(dz,error,'b*-','Linewidth',1.5); hold on
H(2)=loglog(dz,dz,'k--','Linewidth',1.5);
H(3)=loglog(dz,dz.^2,'k-','Linewidth',1.5);

plot_style
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0.01 0.2 0.00001 1])

xlabel('Grid Spacing, $\delta \!z$ (m)','Interpreter','LaTeX')
ylabel('Friction Factor Error')
h = legend(H,'FDS','$O(\delta \!z)$','$O(\delta \!z^2)$','Location','Southeast');
set(h,'Interpreter','LaTeX')
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_5_Verification_Guide/FIGURES/poiseuille_convergence')
