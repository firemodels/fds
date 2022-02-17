% McDermott
% 3-8-13
% synthetic_eddy_method.m

close all
clear all

plot_style

datadir='../../Verification/Turbulence/';
plotdir='../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';

error_tolerance = 0.01;

% Flat profile
% ------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

if ~exist([datadir,'sem_flat_leddy_p2_line.csv'])
   display(['Error: File ',[datadir,'sem_flat_leddy_p2_line.csv'],' does not exist. Skipping case.'])
   return
end

M = importdata([datadir,'sem_flat_leddy_p2_line.csv'],',',2);

k = find(strcmp(M.colheaders,'umean'));
z = M.data(:,k-1);
umean = M.data(:,k);

u0=1;
uprof = u0*ones(length(z),1);
H(1)=plot(uprof,z,'k-'); hold on
H(2)=plot(1.1*uprof,z,'k--');
plot(0.9*uprof,z,'k--')

H(3)=plot(umean,z,'b>-');
axis([.4 1.2 0 1])
xlabel('{\it u} (m/s)','FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\it z} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','northwest');
set(h,'Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Key_Font_Size)

Git_Filename = [datadir,'sem_flat_leddy_p2_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_flat_leddy_p2'])

% compute error
umean_error = norm(umean-uprof)/u0/length(umean);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_flat_leddy_p2.fds umean_error = ',umean_error])
end

urms_error = norm(umean+urms-1.1*uprof)/u0/length(umean);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_flat_leddy_p2.fds urms_error = ',urms_error])
end


% Parabolic profile
% -----------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

if ~exist([datadir,'sem_par_leddy_p2_line.csv'])
   display(['Error: File ',[datadir,'sem_par_leddy_p2_line.csv'],' does not exist. Skipping case.'])
   return
end

M = importdata([datadir,'sem_par_leddy_p2_line.csv'],',',2);

k = find(strcmp(M.colheaders,'umean'));
z = M.data(:,k-1);
umean = M.data(:,k);

umax=1;
h=.5;
a=umax/h^2;
uprof = umax - a*(z-h).^2;
H(1)=plot(uprof,z,'k-'); hold on
H(2)=plot(1.1*uprof,z,'k--');
plot(0.9*uprof,z,'k--')

H(3)=plot(umean,z,'b>-');
axis([0 1.2 0 1])
xlabel('{\it u} (m/s)','FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\itz} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','west');
set(h,'Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

Git_Filename = [datadir,'sem_par_leddy_p2_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_par_leddy_p2'])

% compute error
umean_error = norm(umean-uprof)/umax/length(umean);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_par_leddy_p2.fds umean_error = ',num2str(umean_error)])
end

urms_error = norm(umean+urms-1.1*uprof)/umax/length(umean);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_par_leddy_p2.fds urms_error = ',num2str(urms_error)])
end

% Atmospheric profile
% -------------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

if ~exist([datadir,'sem_atm_leddy_p2_line.csv'])
   display(['Error: File ',[datadir,'sem_atm_leddy_p2_line.csv'],' does not exist. Skipping case.'])
   return
end

M = importdata([datadir,'sem_atm_leddy_p2_line.csv'],',',2);

k = find(strcmp(M.colheaders,'umean'));
z = M.data(:,k-1);
umean = M.data(:,k);

u0=1;
z0=0.5;
p=0.3;
uprof = u0*(z/z0).^p;
H(1)=plot(uprof,z,'k-'); hold on
H(2)=plot(1.1*uprof,z,'k--');
plot(0.9*uprof,z,'k--')

H(3)=plot(umean,z,'b>-');
axis([0 1.4 0 1])
xlabel('{\it u} (m/s)','FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\it z} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','northwest');
set(h,'Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

Git_Filename = [datadir,'sem_atm_leddy_p2_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_atm_leddy_p2'])

% compute error
umean_error = norm(umean-uprof)/u0/length(umean);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_atm_leddy_p2.fds umean_error = ',num2str(umean_error)])
end

urms_error = norm(umean+urms-1.1*uprof)/u0/length(umean);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_atm_leddy_p2.fds urms_error = ',num2str(urms_error)])
end

% RAMP profile
% -------------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

if ~exist([datadir,'sem_ramp_leddy_p2_line.csv'])
   display(['Error: File ',[datadir,'sem_ramp_leddy_p2_line.csv'],' does not exist. Skipping case.'])
   return
end

M = importdata([datadir,'sem_ramp_leddy_p2_line.csv'],',',2);

k = find(strcmp(M.colheaders,'umean'));
z = M.data(:,k-1);
umean = M.data(:,k);

u0=1;
for i=1:length(z)
    if z(i)<.5
        uprof(i) = z(i)*2*u0;
    else
        uprof(i) = (1-z(i))*2*u0;
    end
end

H(1)=plot(uprof,z,'k-'); hold on
H(2)=plot(1.1*uprof,z,'k--');
plot(0.9*uprof,z,'k--')

H(3)=plot(umean,z,'b>-');
axis([0 1.2 0 1])
xlabel('{\it u} (m/s)','FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\it z} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','northeast');
set(h,'Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

Git_Filename = [datadir,'sem_ramp_leddy_p2_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_ramp_leddy_p2'])

% compute error
umean_error = norm(umean-uprof)/u0/length(umean);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_ramp_leddy_p2.fds umean_error = ',num2str(umean_error)])
end

urms_error = norm(umean+urms-1.1*uprof)/u0/length(umean);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_ramp_leddy_p2.fds urms_error = ',num2str(urms_error)])
end

% Monin-Obukhov profile at OPEN inflow boundary VELOCITY
% ------------------------------------------------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

if ~exist([datadir,'sem_open_wind_line.csv'])
   display(['Error: File ',[datadir,'sem_open_wind_line.csv'],' does not exist. Skipping case.'])
   return
end

% expected values
E = importdata([datadir,'sem_open_wind_MO_profile.csv'],',',1);
z_exp = E.data(:,find(strcmp(E.colheaders,'z')));
u_exp = E.data(:,find(strcmp(E.colheaders,'u')));
T_exp = E.data(:,find(strcmp(E.colheaders,'T')));

M = importdata([datadir,'sem_open_wind_line.csv'],',',2);

z_fds = M.data(:,find(strcmp(M.colheaders,'UMEAN-z')));
u_fds = M.data(:,find(strcmp(M.colheaders,'UMEAN')));
u_fds_rms = M.data(:,find(strcmp(M.colheaders,'URMS')));
T_fds = M.data(:,find(strcmp(M.colheaders,'TMEAN')));
T_fds_rms = M.data(:,find(strcmp(M.colheaders,'TRMS')));

I = 0.1; % turbulence intensity (I), VEL_RMS=1 m/s, U_REF = 10 m/s from input file
H(1)=plot(u_exp,z_exp,'k>'); hold on
H(2)=plot((1+I)*u_exp,z_exp,'k:');
plot((1-I)*u_exp,z_exp,'k:')

H(3)=plot(u_fds,z_fds,'b-');
H(4)=plot(u_fds+u_fds_rms,z_fds,'b--');
plot(u_fds-u_fds_rms,z_fds,'b--')
axis([0 12 0 10])
xlabel('{\it u} (m/s)','FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\it z} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

h = legend(H,'Monin-Obukhov profile','Prescribed rms','FDS mean','FDS rms','location','northwest');
set(h,'Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

Git_Filename = [datadir,'sem_open_wind_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_open_wind_u_prof'])

% compute error
u_ref = 10;
kk = find(z_exp<max(z_fds)&z_exp>min(z_fds));
u_fds_int = interp1(z_fds,u_fds,z_exp(kk));
umean_error = norm(u_exp(kk)-u_fds_int)/u_ref/length(u_fds_int);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_open_wind.fds umean_error = ',num2str(umean_error)])
end

u_fds_rms_int = interp1(z_fds,u_fds_rms,z_exp(kk));
urms_error = norm(u_fds_int+u_fds_rms_int-(1+I)*u_exp(kk))/u_ref/length(u_fds_rms_int);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_open_wind.fds urms_error = ',num2str(urms_error)])
end

% Monin-Obukhov profile at OPEN inflow boundary TEMPERATURE
% ---------------------------------------------------------
figure
set(gca,'Units',Plot_Units)
set(gca,'Position',[Plot_X Plot_Y Plot_Width Plot_Height])

H(1)=plot(T_exp,z_exp,'ko'); hold on
H(2)=plot(T_fds,z_fds,'r-');
xlabel('{\it T} (C)','FontName',Font_Name,'FontSize',Label_Font_Size)
ylabel('{\it z} (m)','FontName',Font_Name,'FontSize',Label_Font_Size)

set(gca,'FontName',Font_Name)
set(gca,'FontSize',Label_Font_Size)

h = legend(H(1:2),'Monin-Obukhov profile','FDS mean','location','northeast');
set(h,'Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

Git_Filename = [datadir,'sem_open_wind_git.txt'];
addverstr(gca,Git_Filename,'linear')

set(gcf,'Visible',Figure_Visibility);
set(gcf,'Units',Paper_Units);
set(gcf,'PaperUnits',Paper_Units);
set(gcf,'PaperSize',[Paper_Width Paper_Height]);
set(gcf,'Position',[0 0 Paper_Width Paper_Height]);

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_open_wind_T_prof'])

% compute error
T_ref = 20+273;
kk = find(z_exp<max(z_fds)&z_exp>min(z_fds));
T_fds_int = interp1(z_fds,T_fds,z_exp(kk));
Tmean_error = norm(T_exp(kk)-T_fds_int)/T_ref/length(T_fds_int);
if Tmean_error>error_tolerance
    display(['Matlab Warning: sem_open_wind.fds Tmean_error = ',num2str(Tmean_error)])
end






