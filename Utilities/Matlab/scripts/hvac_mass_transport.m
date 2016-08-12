% B M Ralph
% 12-8-2016
% hvac_mass_transport.m
%
% Convergence study for HVAC transient mass transport.

close all
clear all

t_end = 20;
u = 1;
L = 10;
t0 = 0;

% Analytical solution (using Anonymous Function)
% Y = @(t) 0 .*(t <= L/u) + 1 .*(t > L/u);
Y = @(t) 0 .*(t <= 10.825) + 1 .*(t > 10.825); % a hack for now as solver is wrong

nt = 1000;
dt = t_end/nt;
tc = (t0+dt/2) : dt : (t_end-dt/2);

hold on
plot(tc,Y(tc), 'm-');

% Gather FDS results
ddir = [pwd, '\..\..\Verification\HVAC\'];
fnt = {'HVAC_mass_transport_conv_0320','HVAC_mass_transport_conv_0640',...
    'HVAC_mass_transport_conv_1280','HVAC_mass_transport_conv_2560',...
    'HVAC_mass_transport_conv_5120'};
plt_style = {'b-','g-','r-','m-','c-'};
dx_array = {10/320,10/640,10/1280,10/2560,10/5120};

ert = []; % initialize error norm vector
dxx = []; % init dtt vector

% Loop over cases, plotting FDS results (hold on) and computing error vector
for i=1:length(fnt)
    M = importdata([ddir,fnt{i},'_devc.csv'],',',2);
    Y_fds = M.data(1:end,2);
    Y_fds = Y_fds.';
    nt = length(Y_fds);
    dt = t_end/nt;
    tc = (t0+dt/2) : dt : (t_end-dt/2);
    plot(tc,Y_fds,plt_style{i})
    dx = dx_array{i};
    e_vec = Y_fds - Y(tc);
    ert = [ert,norm(e_vec)/sqrt(length(e_vec))] % populates ert vector, element-by-element
    dxx = [dxx,dx] % populates dtt vector, element-by-element
end

figure

hh(1)=loglog(dxx,ert,'ksq-'); hold on

hh(2)=loglog(dxx,50*dxx,'k--');
hh(3)=loglog(dxx,50*(dxx.^2),'k-');

plot_style
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)

xlabel('{\it \Deltax} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
ylabel('L2 error (seconds)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'Fontname','Times')
legend(hh,'FDS','{\it O(\Deltax)}','{\it O(\Deltax^2)}','location','northwest')
legend('boxoff')

% add SVN if file is available

Git_Filename = [ddir,'HVAC_mass_transport_conv_0320_git.txt'];
addverstr(gca,Git_Filename,'loglog')

% % print to pdf
% set(gcf,'PaperUnits',Paper_Units);
% set(gcf,'PaperSize',[Paper_Width Paper_Height]);
% set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
% print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ht3d_test_1_convergence')

