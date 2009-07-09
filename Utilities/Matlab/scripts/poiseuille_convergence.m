% McDermott
% 5-14-2009
% poiseuille_convergence.m

close all
clear all

paper_width  = 6.0; % inches
paper_height = 4.5; % inches

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

set(gca,'Units','inches')
set(gca,'FontName','Times')
set(gca,'FontSize',14)
set(gca,'Position',[1,0.75,4.5,3.45])

xlabel('$\delta \!z$','Interpreter','LaTeX')
ylabel('friction factor error')
h = legend(H,'FDS','$O(\delta \!z)$','$O(\delta \!z^2)$','Location','Southeast');
set(h,'Interpreter','LaTeX')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[paper_width paper_height]);
set(gcf,'PaperPosition',[0 0 paper_width paper_height]);
print(gcf,'-dpdf','../../../Manuals/FDS_5_Verification_Guide/FIGURES/poiseuille_convergence')
