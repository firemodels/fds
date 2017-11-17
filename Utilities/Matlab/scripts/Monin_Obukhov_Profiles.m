% McGrattan
% 08-17-2017
% Monin_Obukhov_Profiles.m
%
% Atmospheric Boundary Layer profiles, based on M-O theory as described in 
% "Falcon Series Data Report", GRI-89/0138, June 1990.
% These plots are used as illustrations in the FDS User's Guide.

clear all
close all

plot_style

rho_0 = 1.2;
g = 9.81;
cp = 1005;
kappa = 0.41;
p_0 = 100000;

L = -100.;
z_0 = 0.01;
z_r = 2.;
u_r = 5.;
p_r = p_0-rho_0*g*(z_r-z_0);
T_r = 20+273;
theta_0 = (p_0/p_r)^0.285*(T_r+(g/cp)*(z_0-z_r));
u_star = kappa*u_r/log(z_r/z_0);
theta_star = u_star^2*theta_0/(g*kappa*L);

z = [z_0 10*z_0 1 2 3 4 5 6 7 8 9 10 15 20 25 30 50 100];

figure(1)
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

u = (u_star/kappa)*(log(z/z_0) - psi_m);
theta = theta_0 + (theta_star/kappa)*(log(z/z_0) - psi_h);
T = theta*(p_0/p_r)^-0.285 - (g/cp)*(z-z_r);
delta_T = T(1)-T_r;
T = T - delta_T;

plot(u,z,'ko-'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 10 0 30])
text(.05,.90,'$u_{\rm ref}(2 \; \hbox{m})=5$ m/s','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.80,['$L=' num2str(L,'%4.0f\n') '$ m'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.70,['$z_0=' num2str(z_0,'%5.3f\n') '$ m'],'FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
text(.05,.60,'$T_0=20 \;^\circ$C','FontName',Font_Name,'FontSize',Label_Font_Size,'Interpreter','LaTeX','Units','normalized')
xlabel('Velocity (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['vel_' 'L=' num2str(L,'%4.0f\n')])

figure(2)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(T-273,z,'ko-'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([15 25 0 30])
xlabel('Temperature (°C)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',['tmp_' 'L=' num2str(L,'%4.0f\n')])

% Print out the u and T profiles

fid = fopen('Monin_Obukhov_Profile.csv','wt','n');
fprintf(fid,'%s\n','z,u,T');
for j=1:numel(z)
   fprintf(fid,'%5.1f, %5.2f, %5.1f\n',z(j),u(j),T(j)-273);
end
fclose(fid);

