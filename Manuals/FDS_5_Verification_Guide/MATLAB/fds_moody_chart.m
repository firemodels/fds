% McDermott
% 2-24-09
% fds_moody_chart.m

close all
clear all

% % test Colebrook formula
% Re = 1e5;
% RR = 0;
% ff = .1;
% [f_cor,error,iter] = colebrook(Re,RR,ff,1e-3) % correct value is f=0.018

n = 100;
Re = logspace(3.3,8,n);
RR = [0,1e-5,1e-4,1e-3,1e-2,1e-1];
tol = 1e-3;

axes('Fontsize',12)

for i=1:length(RR)

    for j=1:n
        % initial guess for f
        if j>1
            ff = f(j-1);
        else
            ff = .1;
        end
        [f(j),error,iter] = colebrook(Re(j),RR(i),ff,tol);
        M(j,1) = Re(j);
        M(j,i+1) = f(j);
    end
    
    loglog(Re,f,'b-'); hold on
end

Re_DNS = logspace(2,3.3);
f_DNS = (24./Re_DNS);
loglog(Re_DNS,f_DNS,'b-')
axis([1e2 1e8 .005 .2]) % based on MYO
xlabel('Re')
ylabel('friction factor')

repository = '../../../Verification/Flowfields/';

% gather FDS results (laminar)
L = 1;
dpdx = -1;
[f_fds(1),Re_fds(1)] = friction_factor_calc(dpdx,L,[repository,'poiseuille_N64_mu025_devc.csv']);
[f_fds(2),Re_fds(2)] = friction_factor_calc(dpdx,L,[repository,'poiseuille_N64_mu0125_devc.csv']);

loglog(Re_fds,f_fds,'b*')

% gather FDS results (turbulent)
dpdx = -.01;
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-0p01_N8_devc.csv']);loglog(Re,f,'bsq')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-0p01_N16_devc.csv']);loglog(Re,f,'b^')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-0p01_N32_devc.csv']);loglog(Re,f,'bo')

dpdx = -1;
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-1_N8_devc.csv']);loglog(Re,f,'bsq')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-1_N16_devc.csv']);loglog(Re,f,'b^')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-1_N32_devc.csv']);loglog(Re,f,'bo')

dpdx = -100;
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-100_N8_devc.csv']);H(1)=loglog(Re,f,'bsq');
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-100_N16_devc.csv']);H(2)=loglog(Re,f,'b^');
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-100_N32_devc.csv']);H(3)=loglog(Re,f,'bo');

%aa = get(gca,'Position') % 0.1300    0.1381    0.7750    0.7869
set(gca,'Position',[0.1300    0.1381    0.75    0.7869])
text(1.3e8,2e-1,'\epsilon/D','Fontsize',12)
text(1.3e8,.82e-2,num2str(RR(2)),'Fontsize',12)
text(1.3e8,1.2e-2,num2str(RR(3)),'Fontsize',12)
text(1.3e8,1.95e-2,num2str(RR(4)),'Fontsize',12)
text(1.3e8,3.85e-2,num2str(RR(5)),'Fontsize',12)
text(1.3e8,1.02e-1,num2str(RR(6)),'Fontsize',12)
legend(H,'N=8','N=16','N=32','Location','Southwest')
set(gcf,'PaperUnits','inches');
set(gcf,'PaperSize',[6 4.5]);
set(gcf,'PaperPositionMode','auto');
print(gcf,'-dpdf','../FIGURES/fds_moody_chart')

close all
clear all


