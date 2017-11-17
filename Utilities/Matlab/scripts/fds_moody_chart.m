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
RR = [0,1e-4,1e-3,1e-2,1e-1];
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

    if i==1; loglog(Re,f,'k-'); hold on; end
    if i==2; loglog(Re,f,'b-'); end
    if i==3; loglog(Re,f,'r-'); end
    if i==4; loglog(Re,f,'g-'); end
    if i==5; loglog(Re,f,'m-'); end
end

Re_DNS = logspace(2,3.3);
f_DNS = (24./Re_DNS);
loglog(Re_DNS,f_DNS,'k-')

plot_style
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

axis([1e2 1e8 .005 .2]) % based on MYO
xlabel('Re_{\it H}','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('\it f','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size,'Rotation',0.0)

outdir = '../../../out/Moody_Chart/FDS_Output_Files/';

% gather FDS results (laminar)
L = 1;
dpdx = -1;
[f_fds(1),Re_fds(1)] = friction_factor_calc(dpdx,L,[outdir,'poiseuille_N64_mu025_devc.csv']);
[f_fds(2),Re_fds(2)] = friction_factor_calc(dpdx,L,[outdir,'poiseuille_N64_mu0125_devc.csv']);

loglog(Re_fds,f_fds,'k*')

% gather FDS results (turbulent)
mu = 1.84e-5;

dpdx = -.01;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-0p01_N8_devc.csv'],mu);loglog(Re,f,'ksq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-0p01_N16_devc.csv'],mu);loglog(Re,f,'k^')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-0p01_N32_devc.csv'],mu);loglog(Re,f,'ko')

dpdx = -1;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-1_N8_devc.csv'],mu);loglog(Re,f,'ksq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-1_N16_devc.csv'],mu);loglog(Re,f,'k^')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-1_N32_devc.csv'],mu);loglog(Re,f,'ko')

dpdx = -100;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-100_N8_devc.csv'],mu);H(1)=loglog(Re,f,'ksq');
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-100_N16_devc.csv'],mu);H(2)=loglog(Re,f,'k^');
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'moody_dpdx=-100_N32_devc.csv'],mu);H(3)=loglog(Re,f,'ko');

dpdx = -.0001;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p0001_dpdx=-p0001_N8_devc.csv'],mu);loglog(Re,f,'bsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p001_dpdx=-p0001_N8_devc.csv'],mu);loglog(Re,f,'rsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p01_dpdx=-p0001_N8_devc.csv'],mu);loglog(Re,f,'gsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p1_dpdx=-p0001_N8_devc.csv'],mu);loglog(Re,f,'msq')

dpdx = -.01;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p0001_dpdx=-p01_N8_devc.csv'],mu);loglog(Re,f,'bsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p001_dpdx=-p01_N8_devc.csv'],mu);loglog(Re,f,'rsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p01_dpdx=-p01_N8_devc.csv'],mu);loglog(Re,f,'gsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p1_dpdx=-p01_N8_devc.csv'],mu);loglog(Re,f,'msq')

dpdx = -1;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p0001_dpdx=-1_N8_devc.csv'],mu);loglog(Re,f,'bsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p001_dpdx=-1_N8_devc.csv'],mu);loglog(Re,f,'rsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p01_dpdx=-1_N8_devc.csv'],mu);loglog(Re,f,'gsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p1_dpdx=-1_N8_devc.csv'],mu);loglog(Re,f,'msq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p0001_dpdx=-1_N16_devc.csv'],mu);loglog(Re,f,'b^')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p001_dpdx=-1_N16_devc.csv'],mu);loglog(Re,f,'r^')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p01_dpdx=-1_N16_devc.csv'],mu);loglog(Re,f,'g^')
[f,Re] = friction_factor_calc(dpdx,2*L,[outdir,'z0=p02_dpdx=-1_N16_devc.csv'],mu);H(4)=loglog(Re,f,'k>');loglog(Re,f,'g>');

dpdx = -100;
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p0001_dpdx=-100_N8_devc.csv'],mu);loglog(Re,f,'bsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p001_dpdx=-100_N8_devc.csv'],mu);loglog(Re,f,'rsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p01_dpdx=-100_N8_devc.csv'],mu);loglog(Re,f,'gsq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p1_dpdx=-100_N8_devc.csv'],mu);loglog(Re,f,'msq')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p0001_dpdx=-100_N16_devc.csv'],mu);loglog(Re,f,'b^')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p001_dpdx=-100_N16_devc.csv'],mu);loglog(Re,f,'r^')
[f,Re] = friction_factor_calc(dpdx,L,[outdir,'z0=p01_dpdx=-100_N16_devc.csv'],mu);loglog(Re,f,'g^')
[f,Re] = friction_factor_calc(dpdx,2*L,[outdir,'z0=p02_dpdx=-100_N16_devc.csv'],mu);loglog(Re,f,'g>')

text(1.3e8,2e-1,'{\it s/H}','Interpreter',Font_Interpreter,'FontSize',Key_Font_Size,'Fontname',Font_Name)
text(1.3e8,1.2e-2,num2str(RR(2)),'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size,'Fontname',Font_Name)
text(1.3e8,1.95e-2,num2str(RR(3)),'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size,'Fontname',Font_Name)
text(1.3e8,3.85e-2,num2str(RR(4)),'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size,'Fontname',Font_Name)
text(1.3e8,1.02e-1,num2str(RR(5)),'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size,'Fontname',Font_Name)
h = legend(H,'\it N_z=8, H=1','\it N_z=16, H=1','\it N_z=32, H=1','\it N_z=16, H=2','Location','Southwest');
set(h,'Interpreter',Font_Interpreter,'FontSize',Key_Font_Size)

% add VerStr if file is available

Git_Filename = [outdir,'moody_dpdx=-0p01_N8_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[1.1*Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 1.1*Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/fds_moody_chart')


