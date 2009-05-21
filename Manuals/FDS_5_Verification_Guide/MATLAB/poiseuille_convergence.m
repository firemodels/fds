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
axes('Fontsize',12)
dz = L./N;
error = abs(f-24./Re);
H(1)=loglog(dz,error,'b*-','Linewidth',1.5); hold on
H(2)=loglog(dz,dz,'k--','Linewidth',1.5);
H(3)=loglog(dz,dz.^2,'k-','Linewidth',1.5);
xlabel('dz')
ylabel('friction factor error')
legend(H,'FDS','O(dz)','O(dz^2)','Location','Southeast')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[6 4.5]);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpdf','../FIGURES/poiseuille_convergence')

close all
clear all