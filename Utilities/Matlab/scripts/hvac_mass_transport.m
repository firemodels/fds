% B M Ralph
% 12-8-2016
% hvac_mass_transport.m
%
% Convergence study for HVAC transient mass transport (mass fraction at
% downstream duct node).

close all
clear all

% Input parameters
t0 = 0;
t_end = 2;
u = 1;
L = 1;

% Analytical solution (using Anonymous Function)
Y = @(t) 0 .*(t <= L/u) + 1 .*(t > L/u);

% Create mass fraction plot and plot analytical solution
nt = 1000;
dt = t_end/nt;
tc = t0 : dt : t_end;
figure
plot(tc,Y(tc), 'k--');
hold on

% Gather FDS results
ddir = [pwd, '\..\..\Verification\HVAC\'];
fnt = {'HVAC_mass_transport_conv_0020','HVAC_mass_transport_conv_0040',...
    'HVAC_mass_transport_conv_0080','HVAC_mass_transport_conv_0160',...
    'HVAC_mass_transport_conv_0320'};
plt_style = {'b-','g-','r-','m-','c-'};
nc_array = {20,40,80,160,320};
dx_array = {1/20,1/40,1/80,1/160,1/320};

L1e = []; % initialize L1 norm error vector
L2e = []; % initialize L2 norm error vector
dxx = []; % init dxx vector

% Loop over cases, plotting FDS results (hold on) and computing L2 error vector
for i=1:length(fnt)
    M = importdata([ddir,fnt{i},'_devc.csv'],',',2);
    Y_fds = M.data(1:end,2); % FDS species data
    t_fds = M.data(1:end,1); % FDS time
    plot(t_fds,Y_fds,plt_style{i})
    L1e = [L1e,1/nc_array{i}*sum(abs(Y(t_fds)-Y_fds))]; % populates L1 norm error vector, element-by-element
    L2e = [L2e,sqrt(1/nc_array{i}*sum((Y(t_fds)-Y_fds).^2))]; % populates L2 error norm vector, element-by-element
    dxx = [dxx,dx_array{i}]; % populates dxx vector, element-by-element
end

% Mass fraction plot settings
plot_style
set(gca,'Units',Plot_Units)
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

xlabel('Time (s)')
ylabel('Mass fraction (kg/kg)')
legend({'Exact solution','FDS N\_CELLS = 20','FDS N\_CELLS = 40','FDS N\_CELLS = 80',...
    'FDS N\_CELLS = 160','FDS N\_CELLS = 320'},'FontSize',Key_Font_Size,'location','southeast')
legend('boxoff')

Git_Filename = [ddir,'HVAC_mass_transport_conv_0320_git.txt'];
addverstr(gca,Git_Filename,'linear')

% % print to pdf
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/hvac_mass_transport_convergence_1')

% Plot L1 norm error convergence
figure
hh(1)=loglog(dxx,L1e,'ksq-'); 
hold on
hh(2)=loglog(dxx,(dxx.^0.5)/10,'k--');

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Delta x} (m)')
ylabel('L1 error (kg/kg)')
legend(hh,'FDS','{\it O(\Delta x^{1/2})}','location','southeast')
legend('boxoff')

Git_Filename = [ddir,'HVAC_mass_transport_conv_0320_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% % print to pdf
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/hvac_mass_transport_convergence_2')

% Plot L2 error convergence
figure
hh(1)=loglog(dxx,L2e,'ksq-'); 
hold on
hh(2)=loglog(dxx,dxx.^(1/4)/10,'k--');

set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)')
ylabel('L2 error (kg/kg)')
legend(hh,'FDS','{\it O(\Delta x^{1/4})}','location','southeast')
legend('boxoff')

Git_Filename = [ddir,'HVAC_mass_transport_conv_0320_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% % print to pdf
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/hvac_mass_transport_convergence_3')

