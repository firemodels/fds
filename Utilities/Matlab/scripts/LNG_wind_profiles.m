% McGrattan
% 12-14-2016
% wind_profiles.m
%
% Atmospheric Boundary Layer profiles, based on M-O theory as described in
% "Falcon Series Data Report", GRI-89/0138, June 1990.

clear all
close all

outdir = '../../../out/LNG_Dispersion/FDS_Output_Files/';
expdir = '../../../exp/LNG_Dispersion/';
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/LNG_Dispersion/';

plot_style

test = {'Burro3','Burro7','Burro8','Burro9',...
        'Coyote3','Coyote5','Coyote6',...
        'Falcon1','Falcon3','Falcon4',...
        'MaplinSands27','MaplinSands34','MaplinSands35'};

g = 9.81;
cp = 1005;
kappa = 0.41;
p_0 = 1000;

z = [0.001 0.1 0.5 1 1.9 3 4 5 6 7 8 9 10 15 20 25 30];

u_star =     [ 0.255  0.372  0.074  0.252  0.310  0.480  0.210  0.0605  0.3053 0.3694  0.19    0.280   0.315];
theta_star = [-0.532 -0.097  0.026 -0.035 -0.890 -0.520  0.039  0.0577 -0.0175 0.1521 -0.180  -0.0745 -0.0879];
L =          [-9.49  -110.8  16.2  -142.  -8.56  -33.2   82.5   4.963  -422.2  69.38  -14.4   -75.5   -81.15];
p_r =        [ 948    940    941    940.   936    939    942    908.9   900.8  906.3   1000.   1000.   1000.];
z_0 =        [ 0.0002 0.0002 0.0002 0.0002 0.0002 0.0002 0.0002 0.008   0.008  0.008   0.00002 0.00002 0.00002];
z_r =   [1 1 1 1 0.5 0.5 0.5 1 1 1 1.9 1.9 1.9];
i_z_r = [4 4 4 4 3   3   3   4 4 4 5   5   5];

for i=1:13

figure(2*i-1)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

M = importdata([expdir,test{i},'_profile.csv'],',',2);
S = importdata([outdir,test{i},'_line.csv'],',',2);

if L(i)<0;
x = (1-16*z/L(i)).^0.25;
psi_m = 2*log((1+x)/2) + log((1+x.^2)/2) - 2*atan(x) + pi/2;
psi_h = 2*log((1+x.^2)/2);
else
psi_m = -5*z/L(i);
psi_h = psi_m;
end

u = (u_star(i)/kappa)*(log(z/z_0(i)) - psi_m);
T_a = M.data(:,3)+273.15;
theta_r = (p_0/p_r(i))^0.285*(T_a+(g/cp)*(M.data(:,1)-z_r(i)));
theta_0 = theta_r(1) - (theta_star(i)/kappa)*(log(z_r(i)/z_0(i)) - psi_h(i_z_r(i)));
theta = theta_0 + (theta_star(i)/kappa)*(log(z/z_0(i)) - psi_h);
T = theta*(p_0/p_r(i))^-0.285 - (g/cp)*(z-z_r(i));

plot(M.data(:,2),M.data(:,1),'ko'); hold on
plot(u,z,'k-'); hold on
plot(S.data(:,2),S.data(:,1),'k--'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([0 15 0 30])
text(.2,27,[test{i},' Velocity'],'Interpreter',Font_Interpreter,'Fontname',Font_Name,'FontSize',Title_Font_Size)
legend_handle=legend('Measured','M-O Theory','Simulated','Location','SouthEast');
set(legend_handle,'Fontsize',Key_Font_Size);
xlabel('Velocity (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
Git_Filename = [outdir,test{i},'_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,test{i},'_vel'])

figure(2*i)
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

plot(M.data(:,3),M.data(:,1),'ko'); hold on
plot(T-273.15,z,'k-'); hold on
plot(S.data(:,3),S.data(:,1),'k--'); hold on
set(gca,'FontName',Font_Name)
set(gca,'FontSize',Title_Font_Size)
axis([M.data(1,3)-5 M.data(1,3)+5 0 30])
text(M.data(1,3)-4.8,27,[test{i},' Temperature'],'Interpreter',Font_Interpreter,'Fontname',Font_Name,'FontSize',Title_Font_Size)
legend_handle=legend('Measured','M-O Theory','Simulated','Location','SouthEast');
set(legend_handle,'Fontsize',Key_Font_Size);
xlabel('Temperature (\circC)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)
ylabel('Height (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter)

% add Git revision if file is available
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);
print(gcf,'-dpdf',[pltdir,test{i},'_tmp'])

% Print out the FDS RAMP lines for velocity and temperature profiles. This is only needed if the M-O parameters change.

%fid = fopen([test{i},'_ramps.txt'],'wt');
%for j=1:length(z)
%    fprintf(fid,'%s %5.3f %s %5.4f %s\n','&RAMP ID=''u profile'', Z=',z(j),', F=',u(j),' /');
%end
%for j=1:length(z)
%    fprintf(fid,'%s %5.3f %s %5.4f %s\n','&RAMP ID=''T profile'', Z=',z(j),', F=',T(j)/T(i_z_r(i)),' /');
%end
%fclose(fid);

clear M u theta theta_r

end
