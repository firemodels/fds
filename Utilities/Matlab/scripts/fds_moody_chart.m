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

plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
%set(gca,'Position',[0.75,0.75,4.5,3.45])
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gcf,'DefaultLineLineWidth',Line_Width)

axis([1e2 1e8 .005 .2]) % based on MYO
xlabel('Re','Interpreter','LaTeX','FontSize',Label_Font_Size)
ylabel('$f$','Interpreter','LaTeX','FontSize',Label_Font_Size,'Rotation',0.0)

repository = '../../../Validation/Moody_Chart/FDS_Output_Files/';

% gather FDS results (laminar)
L = 1;
dpdx = -1;
[f_fds(1),Re_fds(1)] = friction_factor_calc(dpdx,L,[repository,'poiseuille_N64_mu025_devc.csv']);
[f_fds(2),Re_fds(2)] = friction_factor_calc(dpdx,L,[repository,'poiseuille_N64_mu0125_devc.csv']);

loglog(Re_fds,f_fds,'b*')

% gather FDS results (turbulent)
mu = 1.84e-5;

dpdx = -.01;
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-0p01_N8_devc.csv'],mu);loglog(Re,f,'ksq')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-0p01_N16_devc.csv'],mu);loglog(Re,f,'r^')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-0p01_N32_devc.csv'],mu);loglog(Re,f,'go')

dpdx = -1;
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-1_N8_devc.csv'],mu);loglog(Re,f,'ksq')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-1_N16_devc.csv'],mu);loglog(Re,f,'r^')
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-1_N32_devc.csv'],mu);loglog(Re,f,'go')

dpdx = -100;
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-100_N8_devc.csv'],mu);H(1)=loglog(Re,f,'ksq');
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-100_N16_devc.csv'],mu);H(2)=loglog(Re,f,'r^');
[f,Re] = friction_factor_calc(dpdx,L,[repository,'moody_dpdx=-100_N32_devc.csv'],mu);H(3)=loglog(Re,f,'go');

text(1.3e8,2e-1,'$\varepsilon/D$','Interpreter','LaTeX','FontSize',Key_Font_Size,'Fontname','Times')
text(1.3e8,.82e-2,num2str(RR(2)),'Interpreter','LaTeX','FontSize',Key_Font_Size,'Fontname','Times')
text(1.3e8,1.2e-2,num2str(RR(3)),'Interpreter','LaTeX','FontSize',Key_Font_Size,'Fontname','Times')
text(1.3e8,1.95e-2,num2str(RR(4)),'Interpreter','LaTeX','FontSize',Key_Font_Size,'Fontname','Times')
text(1.3e8,3.85e-2,num2str(RR(5)),'Interpreter','LaTeX','FontSize',Key_Font_Size,'Fontname','Times')
text(1.3e8,1.02e-1,num2str(RR(6)),'Interpreter','LaTeX','FontSize',Key_Font_Size,'Fontname','Times')
h = legend(H,'$N_z=8$','$N_z=16$','$N_z=32$','Location','Southwest');
set(h,'Interpreter','LaTeX')
% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../../Manuals/FDS_5_Verification_Guide/FIGURES/fds_moody_chart')


