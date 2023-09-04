% McDermott
% 15 Nov 2018
% freecon_sphere.m

close all
clear all

plot_style

results_dir = ['../../../out/Convection/'];

% calculations below were used for input file setup

g = 9.80665;
T1 = [303];
T2 = 293;
Tm = 0.5*(T1+T2);
beta = 1./Tm;
MW = 28.85476; % FDS 'LJ AIR'
P0 = 101325;
rho = P0*MW./(8341.5*Tm);
mu = 1.8216e-5;
cp = 1000;
k=0.018216; % for Pr=1 fluid

Pr = cp*mu/k;
nu = mu./rho;
alpha = k./(rho*cp);

setup=1;

if setup

    % see J.P. Holman p. 357 for correlations

    d1 = 0.001;
    Ra1 = (g*beta.*(T1-T2)*d1^3)./(alpha.*nu);
    Nu1 = 2 + .43*Ra1.^(.25);
    Tau1 = 1./Nu1 * d1^2./alpha;

    d2 = .01;
    Ra2 = (g*beta.*(T1-T2)*d2^3)./(alpha.*nu);
    Nu2 = 2 + .43*Ra1.^(.25);
    Tau2 = 1./Nu2 * d2^2./alpha;

    d3 = .1;
    Ra3 = (g*beta.*(T1-T2)*d3^3)./(alpha.*nu);
    Nu3 = 2 + .43*Ra3.^(.25);
    Tau3 = 1./Nu3 * d3^2./alpha;

    d4 = 1;
    Ra4 = (g*beta.*(T1-T2)*d4^3)./(alpha.*nu);
    Nu4 =2 + .5*Ra4.^(.25);
    Tau4 = 1./Nu4 * d4^2./alpha;

    RAYLEIGH_1 = logspace(0,5,100);   % Yuge
    RAYLEIGH_2 = logspace(5,10,100);  % Amato and Tien

    NUSSELT_1 = 2. + 0.43*RAYLEIGH_1.^(.25);
    NUSSELT_2 = 2. + 0.50*RAYLEIGH_2.^(.25);

    % figure(1)
    set(gca,'Units',Plot_Units)
    set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

    marker_handle(1)=loglog(RAYLEIGH_1,NUSSELT_1,'k-'); hold on
    marker_handle(2)=loglog(RAYLEIGH_2,NUSSELT_2,'k--'); hold on
    % loglog(Ra1,Nu1,'k-')
    % loglog(Ra2,Nu2,'k-')
    % loglog(Ra3,Nu3,'k-')
    % loglog(Ra4,Nu4,'k-')

    set(gca,'FontName',Font_Name)
    set(gca,'FontSize',Label_Font_Size)
    axis([1 1e10 1 5e2])
    xlabel('Rayleigh Number','FontSize',Label_Font_Size)
    ylabel('Nusselt Number','FontSize',Label_Font_Size)

    % return
end

% FDS results

casename={...
'free_conv_sphere_1',...
'free_conv_sphere_2',...
'free_conv_sphere_3',...
'free_conv_sphere_4'...
};

delta = [0.001 0.01 0.1 1];
T = [303 303 303 303];
A = 4*pi*(delta/2).^2;

marker_style = {'r^','bsq'};
res = {'8','16'};

for j=1:length(res)
    for i=1:length(delta)

        M = importdata([results_dir,casename{i},'_',res{j},'_devc.csv']);

        % check for steady state
        t = M.data(:,find(strcmp(M.colheaders,'Time')));
        Q = mean(M.data(round(end/2):end,find(strcmp(M.colheaders,'"Q"'))))*1000;
        alpha = k/(rho*cp);
        nu = mu/rho;
        b = 2./(T(i)+T2);
        Ra(i) = (g*b*(T(i)-T2)*delta(i)^3)/(alpha*nu);

        % figure(2)
        % plot(t,Q,'r-'); hold on

        Nu_FDS = (Q/A(i))*(delta(i)/k)/(T2-T(i));

        % figure(1)
        set(gca,'Units',Plot_Units)
        set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])
        marker_handle(2+j)=plot(Ra(i),Nu_FDS,marker_style{j},'MarkerSize',8);
    end
end

lh=legend(marker_handle,'Yuge (1960)','Amato & Tien (1972)','FDS {\itD/\deltax}=8','FDS {\itD/\deltax}=16','Location','Northwest');
set(lh,'FontName',Font_Name,'FontSize',Key_Font_Size)

% title('Free Convection from a Sphere','FontSize',Label_Font_Size)

% add Git if file is available

Git_Filename = [results_dir,'free_conv_sphere_4_16_git.txt'];
addverstr(gca,Git_Filename,'loglog')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperSize',[Paper_Width 1.1*Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/free_conv_sphere');

