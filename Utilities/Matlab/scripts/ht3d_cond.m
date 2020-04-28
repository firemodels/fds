% McDermott
% 6-8-2016
% ht3d_cond.m
%
% Analytical solution to 1D heat equation.
% J. Crank, The Mathematics of Diffusion, 2nd Ed., Oxford, 1975.

close all
clear all

plot_style

t_end = 0.01;
x0 = -.5;
L = 1;
A = 100;
T0 = 20;
lambda = 2*pi/L;
alpha = 1;

T = @(x,t) T0 + A*sin(lambda*(x-x0))*exp(-lambda^2*alpha*t);

nx = 256;
dx = L/nx;
xc = (x0+dx/2):dx:((x0+L)-dx/2);

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

hh(1)=plot(xc,T(xc,0),'k--');
hold on
hh(2)=plot(xc,T(xc,t_end),'k-');

% gather FDS results

ddir = '../../Verification/Heat_Transfer/';
fnx = {'ht3d_nx_10','ht3d_nx_20','ht3d_nx_40','ht3d_nx_80','ht3d_nx_160'};
fny = {'ht3d_ny_10','ht3d_ny_20','ht3d_ny_40','ht3d_ny_80','ht3d_ny_160'};
fnz = {'ht3d_nz_10','ht3d_nz_20','ht3d_nz_40','ht3d_nz_80','ht3d_nz_160'};
plt_style = {'b-o','r-^','g-*','c-sq','m.-'};

erx = []; % initialize error norm vector
dxx = []; % init dxx vector
for i=1:length(fnx)
    M = importdata([ddir,fnx{i},'_devc.csv'],',',2);
    T_fds = M.data(end,2:end);
    nx = length(T_fds);
    dx = L/nx;
    xc = (x0+dx/2):dx:((x0+L)-dx/2);
    plot(xc,T_fds,plt_style{i})
    e_vec = T_fds - T(xc,t_end);
    erx = [erx,norm(e_vec)/sqrt(length(e_vec))];
    dxx = [dxx,dx];
end

ery = []; % initialize error norm vector
dyy = []; % init dxx vector
for i=1:length(fny)
    M = importdata([ddir,fny{i},'_devc.csv'],',',2);
    T_fds = M.data(end,2:end);
    nx = length(T_fds);
    dx = L/nx;
    xc = (x0+dx/2):dx:((x0+L)-dx/2);
    plot(xc,T_fds,plt_style{i})
    e_vec = T_fds - T(xc,t_end);
    ery = [ery,norm(e_vec)/sqrt(length(e_vec))];
    dyy = [dyy,dx];
end

erz = []; % initialize error norm vector
dzz = []; % init dxx vector
for i=1:length(fnz)
    M = importdata([ddir,fnz{i},'_devc.csv'],',',2);
    T_fds = M.data(end,2:end);
    nx = length(T_fds);
    dx = L/nx;
    xc = (x0+dx/2):dx:((x0+L)-dx/2);
    hh(i+2)=plot(xc,T_fds,plt_style{i});
    e_vec = T_fds - T(xc,t_end);
    erz = [erz,norm(e_vec)/sqrt(length(e_vec))];
    dzz = [dzz,dx];
end

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it x} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('{\it T} (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
lh=legend(hh,'Initial Condition','Final Exact','FDS {\itnx}=10','FDS {\itnx}=20','FDS {\itnx}=40','FDS {\itnx}=80','FDS {\itnx}=160','location','northeast');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add version string if file is available

Git_Filename = [ddir,'ht3d_nx_160_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_test_1_profile')

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

clear hh

hh(1)=loglog(dxx,erx,'k+-'); hold on
hh(2)=loglog(dyy,ery,'rsq-');
hh(3)=loglog(dzz,erz,'bo-');

error_tol = 0.005;
if min(erx)>error_tol
    display(['Error: ht3d_nx out of tolerance'])
    erx
end
if min(ery)>error_tol
    display(['Error: ht3d_ny out of tolerance'])
    ery
end
if min(erz)>error_tol
    display(['Error: ht3d_nz out of tolerance'])
    erz
end

hh(4)=loglog(dxx,100*dxx,'k--');
hh(5)=loglog(dxx,500*dxx.^2,'k-');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('L2 error (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
lh=legend(hh,'FDS {\itnx}','FDS {\itny}','FDS {\itnz}','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add version string if file is available

Git_Filename = [ddir,'ht3d_nx_160_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_test_1_convergence')

