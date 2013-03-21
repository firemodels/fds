% McDermott
% 3-8-13
% synthetic_eddy_method.m

close all
clear all

datadir='../../Verification/Turbulence/';
plotdir='../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/';

plot_style

error_tolerance = 0.01;

% Flat profile
% ------------
figure

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
xlabel('$u$ (m/s)')
ylabel('$z$ (m)')

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','northwest');
set(h,'Interpreter',Font_Interpreter)

% add SVN if file is available

SVN_Filename = [datadir,'sem_flat_leddy_p2_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

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
xlabel('$u$ (m/s)')
ylabel('$z$ (m)')

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','west');
set(h,'Interpreter',Font_Interpreter)

% add SVN if file is available

SVN_Filename = [datadir,'sem_par_leddy_p2_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_par_leddy_p2'])

% compute error
umean_error = norm(umean-uprof)/umax/length(umean);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_par_leddy_p2.fds umean_error = ',umean_error])
end

urms_error = norm(umean+urms-1.1*uprof)/umax/length(umean);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_par_leddy_p2.fds urms_error = ',urms_error])
end

% Atmospheric profile
% -------------------
figure

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
xlabel('$u$ (m/s)')
ylabel('$z$ (m)')

k = find(strcmp(M.colheaders,'urms'));
z = M.data(:,k-1);
urms = M.data(:,k);

H(4)=plot(umean+urms,z,'r--');
plot(umean-urms,z,'r--')

h = legend(H,'Prescribed mean','Prescribed rms','FDS mean','FDS rms','location','northwest');
set(h,'Interpreter',Font_Interpreter)

% add SVN if file is available

SVN_Filename = [datadir,'sem_atm_leddy_p2_svn.txt'];
if exist(SVN_Filename,'file')
    SVN = importdata(SVN_Filename);
    x_lim = get(gca,'XLim');
    y_lim = get(gca,'YLim');
    X_SVN_Position = x_lim(1)+SVN_Scale_X*(x_lim(2)-x_lim(1));
    Y_SVN_Position = y_lim(1)+SVN_Scale_Y*(y_lim(2)-y_lim(1));
    text(X_SVN_Position,Y_SVN_Position,['SVN ',num2str(SVN)], ...
        'FontSize',10,'FontName',Font_Name,'Interpreter',Font_Interpreter)
end

% print to pdf
print(gcf,'-dpdf',[plotdir,'sem_atm_leddy_p2'])

% compute error
umean_error = norm(umean-uprof)/u0/length(umean);
if umean_error>error_tolerance
    display(['Matlab Warning: sem_atm_leddy_p2.fds umean_error = ',umean_error])
end

urms_error = norm(umean+urms-1.1*uprof)/u0/length(umean);
if urms_error>error_tolerance
    display(['Matlab Warning: sem_atm_leddy_p2.fds urms_error = ',urms_error])
end







