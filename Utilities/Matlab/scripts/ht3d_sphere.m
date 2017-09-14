% Salah BENKORICHI & Randy McDermott
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

a1   = 0.10334;     % m
a2   = 0.10295;     % m
a3   = 0.11;        % m

n1 = 32;
n2 = 22;
n3 = 12;
r1 = linspace(0,a1,n1);
r2 = linspace(0,a2,n2);
r3 = linspace(0,a3,n3);

t = [10 20 60 120 180]; % seconds

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

for m=1:length(t)

    for x=1:length(r1)

        sum_term = 0;
        for n=1:n1
            sum_term = sum_term + (-1)^n/n^3 * sin(n*pi*r1(x)/a1) * exp(-alpha*t(m)*(n*pi/a1)^2);
        end

        DT1(x) = 20 + g0/(6*k) * (a1^2 - r1(x)^2) + 2*g0*a1^3/(k*pi^3*r1(x)) * sum_term;

    end

    Exact = plot(r1(1:1:end),DT1(1:1:end),'k-x'); hold on

    for y=1:length(r2)

        sum_term = 0;
        for n=1: n2
            sum_term = sum_term + (-1)^n/n^3 * sin(n*pi*r2(y)/a2) * exp(-alpha*t(m)*(n*pi/a2)^2);
        end

        DT2(y) = 20 + g0/(6*k) * (a2^2 - r2(y)^2) + 2*g0*a2^3/(k*pi^3*r2(y)) * sum_term;

    end

    %Exact = plot(r2,DT2,'k'); hold on

    for w=1:length(r3)

        sum_term = 0;
        for n=1:n3
            sum_term = sum_term + (-1)^n/n^3 * sin(n*pi*r3(w)/a3) * exp(-alpha*t(m)*(n*pi/a3)^2);
        end

        DT3(w) = 20 + g0/(6*k) * (a3^2 - r3(w)^2) + 2*g0*a3^3/(k*pi^3*r3(w)) * sum_term;

    end

    %Exact = plot(r3,DT3,'k'); hold on

end

%% gather FDS results

ddir = '../../Verification/Heat_Transfer/';
fnt = {'ht3d_sphere_75'};
fileName = {'ht3d_sphere_25','ht3d_sphere_51','ht3d_sphere_75'};
nc_array = {25,51,75};
dx_array = {0.25/25,0.25/51,0.25/75};

M = importdata([ddir,fnt{1},'_devc.csv'],',',2);
T_fds1 = M.data(2,2:end);
T_fds2 = M.data(3,2:end);
T_fds3 = M.data(7,2:end);
T_fds4 = M.data(13,2:end);
T_fds5 = M.data(end,2:end);
t10 =  plot(r1, T_fds1,'--ob');
t20 =  plot(r1, T_fds2,'--og');
t60 =  plot(r1, T_fds3,'--or');
t120 = plot(r1, T_fds4,'--oc');
t180 = plot(r1, T_fds5,'--om');

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

axis([0 0.105 20 60])
xlabel('Radial Distance (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('Temperature (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
lh=legend([Exact, t10, t20, t60, t120, t180], {'Analytical', 'FDS {\itt}=10 s', 'FDS {\itt}=20 s', 'FDS {\itt}=60 s', 'FDS {\itt}=120 s', 'FDS {\itt}=180 s'},'location','west');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add version string if file is available

Git_Filename = [ddir,'ht3d_sphere_75_git.txt'];
addverstr(gca,Git_Filename,'linear')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_sphere_profile')

% estimating L1 & L2 norm errors

L1e = []; % initialize L1 norm error vector
L2e = []; % initialize L2 norm error vector
dxx = []; % init dxx vector
DT  = {DT3(2:end),DT2(2:end),DT1(2:end)} ;
r   = {r3(2:end),r2(2:end),r1(2:end)};

for i=1:length(fileName)
    M1 = importdata([ddir,fileName{i},'_devc.csv'],',',2);
    T_fds = M1.data(end,2:end); % FDS devices data
    T_fds1=T_fds(2:end);
    L1e = [L1e,1/nc_array{i}*sum(abs(DT{i}-T_fds1))]; % populates L1 norm error vector, element-by-element
    L2e = [L2e,sqrt(1/nc_array{i}*sum((DT{i}-T_fds1).^2))]; % populates L2 error norm vector, element-by-element
    dxx = [dxx,dx_array{i}]; % populates dxx vector, element-by-element
end

% Set the figure

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

hh(1)=loglog(dxx,L1e,'msq-'); hold on
hh(2)=loglog(dxx,5e1*dxx,'k--');
hh(3)=loglog(dxx,5e3*dxx.^2,'k-');

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('L1 error (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
lh=legend(hh,'FDS','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% add version string if file is available

Git_Filename = [ddir,'ht3d_sphere_75_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_sphere_convergence1')

% Plot L2 error convergence

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

clear hh

hh(1)=loglog(dxx,L2e,'rsq-'); hold on
hh(2)=loglog(dxx,2e2*dxx,'k--');
hh(3)=loglog(dxx,2e4*dxx.^2,'k-');
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)')
ylabel('L2 error (\circC)')
lh=legend(hh,'FDS','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

Git_Filename = [ddir,'ht3d_sphere_75_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% print to pdf
set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_sphere_convergence2')