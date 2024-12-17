% McGrattan
% 08-17-2017
% atmospheric_boundary_layer.m
%
% Atmospheric Boundary Layer profiles, based on M-O theory as described in 
% "Falcon Series Data Report", GRI-89/0138, June 1990.
% These plots are used as illustrations in the FDS User's Guide.

clear all
close all

plot_style

basein ='../../Verification/Atmospheric_Effects/atmospheric_boundary_layer';
baseout='../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/atmospheric_boundary_layer';

for i=1:4

datafile=[basein '_' num2str(i,'%1d\n')];
outfile =[baseout '_'  num2str(i,'%1d\n')];

M1 = importdata([datafile,'_devc.csv'],',',2);
M2 = importdata([datafile,'_line.csv'],',',2);

u_r = M1.data(end,2);
T_r = M1.data(end,3) + 273.15;

rho_0 = 1.2;
g = 9.81;
cp = 1005;
kappa = 0.41;
p_0 = 100000;
qdot = {50 -50 25 -25};
z_0 = {0.25 0.25 0.125 0.125};
T_low  = {15 15 15 15};
T_high = {25 25 25 25};
u_high = {20 20  10 15};
fvec = {0.01 0.01 0.002 0.005};
s = {8.15 8.15 4.075 4.075};

theta_0 = T_r;
z_r = 10.;
p_r = p_0-rho_0*g*(z_r-z_0{i});
theta_r = T_r*(p_0/p_r)^0.285;
u_star = kappa*u_r/log(z_r/z_0{i});
L = -u_star^3*theta_0*rho_0*cp/(g*kappa*qdot{i});
theta_star = u_star^2*theta_0/(g*kappa*L);

z = [z_0{i} 10*z_0{i} 1 2 3 4 5 6 7 8 9 10 15 20 25 30 50 100];

f1=figure;
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

if L<0;
x = (1-16*z/L).^0.25;
psi_m = 2*log((1+x)/2) + log((1+x.^2)/2) - 2*atan(x) + pi/2;
psi_h = 2*log((1+x.^2)/2);
else
psi_m = -5*z/L;
psi_h = psi_m;
end

u = (u_star/kappa)*(log(z/z_0{i}) - psi_m);
theta = theta_0 + (theta_star/kappa)*(log(z/z_0{i}) - psi_h);
T = theta.*(p_0./(p_0-rho_0*g*(z-z_0{i}))).^-0.285;

T = T + (theta_0-T(12));

ERROR = abs(u(end)-M2.data(end,2));
if ERROR>2.
   display(['Matlab Warning: atmospheric_boundary_layer Case ',num2str(i),' velocity out of tolerance. ERROR = ',num2str(ERROR),' m/s'])
end

plot(u,z,'ko'); hold on
plot(M2.data(:,2),M2.data(:,1),'k-'); hold off
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 u_high{i} 0 100])
text(.05,.90,['Case ' num2str(i,'%1d\n')],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')

text(.05,.80,['$F=' num2str(fvec{i},'%5.3f\n') '\; \hbox{Pa/m}$'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.70,['$s=' num2str(s{i},'%5.2f\n') '\; \hbox{m}$'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.60,['$\dot{q}_{\rm c}^{\prime \prime}=' num2str(qdot{i}/1000,'%5.3f\n') '\; \hbox{kW/m}^2$'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')

text(.05,.50,['$u(' num2str(z_r,'%2.0f\n') '\; \hbox{m})=' num2str(u_r,'%4.1f\n') '\; \hbox{m/s}$'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.40,['$L=' num2str(L,'%4.0f\n') '$ m'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.30,['$z_0=' num2str(z_0{i},'%5.3f\n') '$ m'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
xlabel('Velocity (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
legend({'M-O Theory','FDS'}, 'Location', 'SouthEast');

git_file = [datafile, '_git.txt'];
addverstr(gca,git_file,'linear')

set(f1,'Visible',Figure_Visibility);
set(f1,'Units',Paper_Units);
set(f1,'PaperUnits',Paper_Units);
set(f1,'PaperSize',[Paper_Width Paper_Height]);
set(f1,'Position',[0 0 Paper_Width Paper_Height]);
print(f1,'-dpdf',[outfile '_vel'])

f2=figure;
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

ERROR = abs(T(end)-273.15-M2.data(end,3));
if ERROR>1.0
   display(['Matlab Warning: atmospheric_boundary_layer Case ',num2str(i),' temperature out of tolerance. ERROR = ',num2str(ERROR),' K'])
end

plot(T-273.15,z,'ko'); hold on
plot(M2.data(:,3),M2.data(:,1),'k-'); hold off
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([T_low{i} T_high{i} 0 100])
text(.05,.90,['Case ' num2str(i,'%1d\n')],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.80,['$T(' num2str(z_r,'%2.0f\n') '\; \hbox{m})=' num2str(T_r-273,'%3.1f\n') '\;^\circ$C'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
xlabel('Temperature (°C)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
legend({'M-O Theory','FDS'}, 'Location', 'SouthWest');

addverstr(gca,git_file,'linear')

set(f2,'Visible',Figure_Visibility);
set(f2,'Units',Paper_Units);
set(f2,'PaperUnits',Paper_Units);
set(f2,'PaperSize',[Paper_Width Paper_Height]);
set(f2,'Position',[0 0 Paper_Width Paper_Height]);
print(f2,'-dpdf',[outfile '_tmp'])

end
