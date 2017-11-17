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






