% McDermott
% 6-8-2016
% crank_1.m
%
% Analytical solution to 1D heat equation.
% J. Crank, The Mathematics of Diffusion, 2nd Ed., Oxford, 1975.

close all
clear all

t_end = 0.01;
x0 = -.5;
L = 1;
A = 100;
T0 = 20;
lamda = 2*pi/L;
alpha = 1;

T = @(x,t) T0 + A*sin(lamda*(x-x0))*exp(-lamda^2*alpha*t);

nx = 256;
dx = L/nx;
xc = (x0+dx/2):dx:((x0+L)-dx/2);

plot(xc,T(xc,0),'k:')
hold on
plot(xc,T(xc,t_end),'k-')

% gather FDS results

ddir = '/Volumes/rmcdermo/GitHub/fds-smv_rmcdermo/Verification/Heat_Transfer/';
fnx = {'ht3d_nx_10','ht3d_nx_20','ht3d_nx_40','ht3d_nx_80','ht3d_nx_160'};
fny = {'ht3d_ny_10','ht3d_ny_20','ht3d_ny_40','ht3d_ny_80','ht3d_ny_160'};
fnz = {'ht3d_nz_10','ht3d_nz_20','ht3d_nz_40','ht3d_nz_80','ht3d_nz_160'};
plt_style = {'b--','r--','g--','c--','m--'};

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
for i=1:length(fny)
    M = importdata([ddir,fnz{i},'_devc.csv'],',',2);
    T_fds = M.data(end,2:end);
    nx = length(T_fds);
    dx = L/nx;
    xc = (x0+dx/2):dx:((x0+L)-dx/2);
    plot(xc,T_fds,plt_style{i})
    e_vec = T_fds - T(xc,t_end);
    erz = [erz,norm(e_vec)/sqrt(length(e_vec))];
    dzz = [dzz,dx];
end

figure

hh(1)=loglog(dxx,erx,'k+-'); hold on
hh(2)=loglog(dyy,ery,'rsq-');
hh(3)=loglog(dzz,erz,'bo-');

hh(4)=loglog(dxx,100*dxx,'k--');
hh(5)=loglog(dxx,500*dxx.^2,'k-');

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('L2 error (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
legend(hh,'FDS nx','FDS ny','FDS nz','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest')
legend('boxoff')

% add SVN if file is available

SVN_Filename = [ddir,'ht3d_nx_160_git.txt'];
addverstr(gca,SVN_Filename,'loglog')

% % print to pdf
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf','../../FDS/Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_test_1_convergence')

