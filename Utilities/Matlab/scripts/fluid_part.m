% McDermott
% 10-12-11
% fluid_part.m

close all
clear all

% set the plot style parameters

ddir='../../Verification/Sprinklers_and_Sprays/';

skip_case = 0;

if ~exist([ddir,'fluid_part_mom_x_devc.csv'])
    display(['Error: File ' [ddir,'fluid_part_mom_x_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([ddir,'fluid_part_mom_y_devc.csv'])
    display(['Error: File ' [ddir,'fluid_part_mom_y_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([ddir,'fluid_part_mom_z_devc.csv'])
    display(['Error: File ' [ddir,'fluid_part_mom_z_devc.csv'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if ~exist([ddir,'fluid_part_mom_x.prt5'])
    display(['Error: File ' [ddir,'fluid_part_mom_x.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([ddir,'fluid_part_mom_y.prt5'])
    display(['Error: File ' [ddir,'fluid_part_mom_y.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([ddir,'fluid_part_mom_z.prt5'])
    display(['Error: File ' [ddir,'fluid_part_mom_z.prt5'] ' does not exist. Skipping case.'])
    skip_case = 1;
end

if skip_case
    return
end

M = importdata([ddir,'fluid_part_mom_x_devc.csv'],',',2);
tx = M.data(:,1);
U = M.data(:,2);
MX = M.data(:,3);

M = importdata([ddir,'fluid_part_mom_y_devc.csv'],',',2);
ty = M.data(:,1);
V = M.data(:,2);
MY = M.data(:,3);

M = importdata([ddir,'fluid_part_mom_z_devc.csv'],',',2);
tz = M.data(:,1);
W = M.data(:,2);
MZ = M.data(:,3);

range=1:10:min([length(tx),length(ty),length(tz)]);

plot_style
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=plot(tx(range),MX(range).*U(range),'bo'); hold on
H(2)=plot(ty(range),MY(range).*V(range),'bv');
H(3)=plot(tz(range),MZ(range).*W(range),'b+');

n = 1000;                 % number of particles
mpv = 10.;                % kg/m^3, mass_per_volume (from FDS input file)
v_xb = 1^3;               % volume of XB region on init line in FDS input file
rho_p = 1000;             % density of water, kg/m^3
d_p = 1000e-6;            % diameter, m
v_p = 4/3*pi*(d_p/2)^3;   % volume of a single droplet, m^3
m_p = rho_p*v_p;          % mass of single droplet, kg
pwt = mpv*v_xb/(n*m_p);   % particle weight factor

[STIME_X, XP_X, YP_X, ZP_X, QP_X] = read_prt5([ddir,'fluid_part_mom_x.prt5'],'real*4');
[STIME_Y, XP_Y, YP_Y, ZP_Y, QP_Y] = read_prt5([ddir,'fluid_part_mom_y.prt5'],'real*4');
[STIME_Z, XP_Z, YP_Z, ZP_Z, QP_Z] = read_prt5([ddir,'fluid_part_mom_z.prt5'],'real*4');

P_X = zeros(1,numel(STIME_X));
for i=1:numel(STIME_X)
    for j=1:n
        P_X(i) = P_X(i) + pwt*m_p*QP_X(i,j,1,1); % momentum of particle at time STIME(i)
    end
end

P_Y = zeros(1,numel(STIME_Y));
for i=1:numel(STIME_Y)
    for j=1:n
        P_Y(i) = P_Y(i) + pwt*m_p*QP_Y(i,j,1,1);
    end
end

P_Z = zeros(1,numel(STIME_Z));
for i=1:numel(STIME_Z)
    for j=1:n
        P_Z(i) = P_Z(i) + pwt*m_p*QP_Z(i,j,1,1);
    end
end

range=1:10:min([length(STIME_X),length(STIME_Y),length(STIME_Z)]);

H(4)=plot(STIME_X(range),P_X(range),'ro');
H(5)=plot(STIME_Y(range),P_Y(range),'rv');
H(6)=plot(STIME_Z(range),P_Z(range),'r+');

P_total_X = P_X + (MX.*U)'; % total momentum
P_total_Y = P_Y + (MY.*V)';
P_total_Z = P_Z + (MZ.*W)';

H(7)=plot(STIME_X(range),P_total_X(range),'go');
H(8)=plot(STIME_Y(range),P_total_Y(range),'gv');
H(9)=plot(STIME_Z(range),P_total_Z(range),'g+');

% analytical solution

rho = 1.1992661;    % fluid density
Cd = 1;             % drag coefficient
A_p = pi*(d_p/2)^2; % particle area

u_p = QP_Z(1,1,1,1);  % initial velocity
U_p = U(1);

M_p = MX(1)/(n*pwt);
alpha = M_p/m_p;


for i=1:(numel(tx)-1)

    u_soln(i) = u_p;
    U_soln(i) = U_p;
    t_soln(i) = tx(i);

    dt = tx(i+1)-tx(i);

    u0 = u_p;
    U0 = U_p;

    beta = 0.5*rho*Cd*A_p*(1/m_p + 1/M_p)*abs(u0-U0);

    u_p = u0/(1+beta*dt) + (u0+alpha*U0)/(1+alpha)*(beta*dt)/(1+beta*dt);

    U_p = U0 + n*pwt*m_p/MX(1)*(u0-u_p);

end

H(10)=plot(t_soln,MX(1)*U_soln,'b-');
H(11)=plot(t_soln,n*pwt*m_p*u_soln,'r-');

axis([min(tx) max(tx) 0 150])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Momentum (kg m/s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

legend_handle = legend(H(1:11),'FDS fluid U','FDS fluid V','FDS fluid W', ...
              'FDS particle U','FDS particle V','FDS particle W', ...
              'FDS total U','FDS total V','FDS total W', ...
              'Analytical fluid','Analytical particle','Location','EastOutside');
set(legend_handle,'FontSize',Key_Font_Size)
set(legend_handle,'Units',Paper_Units)
pos = get(legend_handle,'position');
set(legend_handle,'position',[Paper_Width pos(2:4)])

% add Git revision if file is available
Git_Filename = [ddir,'fluid_part_mom_x_git.txt'];
addverstr(gca,Git_Filename,'linear')

PDF_Paper_Width = 1.5*Paper_Width;

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 PDF_Paper_Width Paper_Height]);

print -dpdf ../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/fluid_part_momentum

% plot velocities

figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=plot(tx(range),U(range),'bo'); hold on
H(2)=plot(ty(range),V(range),'bv');
H(3)=plot(tz(range),W(range),'b+');

U_p = P_X/(n*pwt*m_p);
V_p = P_Y/(n*pwt*m_p);
W_p = P_Z/(n*pwt*m_p);

H(4)=plot(tx(range),U_p(range),'ro');
H(5)=plot(ty(range),V_p(range),'rv');
H(6)=plot(tz(range),W_p(range),'r+');

U_eq = P_total_X/(MX(1)+n*pwt*m_p);
V_eq = P_total_Y/(MY(1)+n*pwt*m_p);
W_eq = P_total_Z/(MZ(1)+n*pwt*m_p);

H(7)=plot(tx,U_eq,'g--');
H(7)=plot(ty,V_eq,'g--');
H(7)=plot(tz,W_eq,'g--');

axis([min(tx) max(tx) 0 10])
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)
xlabel('Time (s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
ylabel('Velocity (m/s)','Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

legend_handle = legend(H(1:7),'FDS fluid U','FDS fluid V','FDS fluid W',...
              'FDS particle U','FDS particle V','FDS particle W',...
              'Equilibrium velocity','Location','EastOutside');
set(legend_handle,'FontSize',Key_Font_Size)
set(legend_handle,'Units',Paper_Units)
pos = get(legend_handle,'position');
set(legend_handle,'position',[Paper_Width pos(2:4)])

% add Git revision if file is available
Git_Filename = [ddir,'fluid_part_mom_x_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[PDF_Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 PDF_Paper_Width Paper_Height]);

print -dpdf ../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/fluid_part_velocity






