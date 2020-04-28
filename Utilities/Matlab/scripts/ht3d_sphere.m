% Salah Benkorichi & Randy McDermott
% 17-7-2017
% ht3d_sphere.m
%
% The solution is from Carslaw and Jaeger, Sec. 9.8, p. 243, Eq. (6)
%
% The notation adopted here follows C. Lautenberger, IAFSS, 2014.

close all
clear all

plot_style

% analytical solution

k     = 1.0;        % W/m/k
rho   = 1000;       % kg/m3
cp    = 1000;       % J/kg/K
g0    = 2e5;        % W/m3
alpha = k/(rho*cp); % m2/s

a1   = 0.10025; % m (this should match the last cell face in ht3d_sphere_102.fds)
a2   = 0.1005;  % m (this should match the last cell face in ht3d_sphere_51.fds)
a3   = 0.105;   % m (this should match the last cell face in ht3d_sphere_25.fds)

n1 = 41;
n2 = 21;
n3 = 11;
r1 = linspace(0.0,0.099,n1); % this should match line DEVC in ht3d_sphere_102.fds
r2 = linspace(0.0,0.098,n2); % this should match line DEVC in ht3d_sphere_51.fds
r3 = linspace(0.0,0.100,n3); % this should match line DEVC in ht3d_sphere_25.fds

t = [10 20 60 120 180]; % seconds

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

for m=1:length(t)

    for ii=1:length(r1)

        sum_term = 0;
        for n=1:n1
            sum_term = sum_term + (-1)^n/n^3 * sin(n*pi*r1(ii)/a1) * exp(-alpha*t(m)*(n*pi/a1)^2);
        end

        DT1(ii) = 20 + g0/(6*k) * (a1^2 - r1(ii)^2) + 2*g0*a1^3/(k*pi^3*r1(ii)) * sum_term;

    end

    % Exact = plot(r1,DT1,'k-x'); hold on

    for jj=1:length(r2)

        sum_term = 0;
        for n=1: n2
            sum_term = sum_term + (-1)^n/n^3 * sin(n*pi*r2(jj)/a2) * exp(-alpha*t(m)*(n*pi/a2)^2);
        end

        DT2(jj) = 20 + g0/(6*k) * (a2^2 - r2(jj)^2) + 2*g0*a2^3/(k*pi^3*r2(jj)) * sum_term;

    end

    Exact = plot(r2,DT2,'k-x'); hold on

    for kk=1:length(r3)

        sum_term = 0;
        for n=1:n3
            sum_term = sum_term + (-1)^n/n^3 * sin(n*pi*r3(kk)/a3) * exp(-alpha*t(m)*(n*pi/a3)^2);
        end

        DT3(kk) = 20 + g0/(6*k) * (a3^2 - r3(kk)^2) + 2*g0*a3^3/(k*pi^3*r3(kk)) * sum_term;

    end

    % Exact = plot(r3,DT3,'k-x'); hold on

end

%% gather FDS results

ddir = '../../Verification/Heat_Transfer/';
fnt = {'ht3d_sphere_51'};
fileName = {'ht3d_sphere_25','ht3d_sphere_51','ht3d_sphere_102'};
nc_array = [25,51,102];
dx_array = 0.25./nc_array;

M = importdata([ddir,fnt{1},'_devc.csv'],',',2);
T_fds1 = M.data(2,2:end);
T_fds2 = M.data(3,2:end);
T_fds3 = M.data(7,2:end);
T_fds4 = M.data(13,2:end);
T_fds5 = M.data(end,2:end);
t10 =  plot(r2,T_fds1,'--ob');
t20 =  plot(r2,T_fds2,'--og');
t60 =  plot(r2,T_fds3,'--or');
t120 = plot(r2,T_fds4,'--oc');
t180 = plot(r2,T_fds5,'--om');

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

axis([0 0.105 20 60])
xlabel('Radial Distance (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('Temperature (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
lh=legend([Exact, t10, t20, t60, t120, t180], {'Analytical', 'FDS {\itt}=10 s', 'FDS {\itt}=20 s', 'FDS {\itt}=60 s', 'FDS {\itt}=120 s', 'FDS {\itt}=180 s'},'location','west');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add version string if file is available

Git_Filename = [ddir,'ht3d_sphere_51_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_sphere_profile')

% estimating L1 & L2 norm errors

Linfe = []; % initialize L1 norm error vector
L2e = []; % initialize L2 norm error vector
dxx = []; % init dxx vector
DT  = {DT3(end),DT2(end),DT1(end)} ;
r   = {r3(end),r2(end),r1(end)};

for i=1:length(fileName)
    M1 = importdata([ddir,fileName{i},'_devc.csv'],',',2);
    T_fds = M1.data(end,end); % FDS devices data
    T_fds1=T_fds(end);
    Linfe = [Linfe,1/(length(T_fds1))*max(abs(DT{i}-T_fds1))]; % populates Linf norm error vector, element-by-element
    L2e = [L2e,sqrt(1/(length(T_fds1))*sum((DT{i}-T_fds1).^2))]; % populates L2 error norm vector, element-by-element
    dxx = [dxx,dx_array(i)]; % populates dxx vector, element-by-element
end

% Set the figure

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

hh(1)=loglog(dxx,Linfe,'msq-'); hold on
hh(2)=loglog(dxx,2e2*dxx,'k--');
hh(3)=loglog(dxx,2e4*dxx.^2,'k-');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
ylabel('{\it L}_{\infty} error (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname',Font_Name)
lh=legend(hh,'FDS','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add version string if file is available

Git_Filename = [ddir,'ht3d_sphere_51_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_sphere_convergence1')

% % Plot L2 error convergence
%%% this is bad case for L2 because the errors are almost nothing on the interior and this skews the results

% figure
% set(gca,'Units',Plot_Units)
% set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

% clear hh

% hh(1)=loglog(dxx,L2e,'rsq-'); hold on
% hh(2)=loglog(dxx,1e2*dxx,'k--');
% hh(3)=loglog(dxx,1e4*dxx.^2,'k-');
% set(gca,'FontName',Font_Name)
% set(gca,'FontSize',Title_Font_Size)

% xlabel('{\it \Deltax} (m)')
% ylabel('L2 error (\circC)')
% lh=legend(hh,'FDS','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest');
% set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% Git_Filename = [ddir,'ht3d_sphere_51_git.txt'];
% addverstr(gca,Git_Filename,'loglog')

% % print to pdf
% set(gcf,'Visible',Figure_Visibility);
% set(gcf,'Units',Paper_Units);
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_sphere_convergence2')
