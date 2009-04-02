% McDermott
% 10-13-2008
% ns2d_convergence.m
%
% This script processes the output for the ns2d analytical solution
% verification test.
%
% The necessary FDS input files in \Verification\NS_Analytical_Solution\ are:
% 
% ns2d_8.fds
% ns2d_16.fds
% ns2d_32.fds
% ns2d_64.fds
% ns2d_8_nupt1.fds
% ns2d_16_nupt1.fds
% ns2d_32_nupt1.fds
% ns2d_64_nupt1.fds


close all
clear all

L = 2*pi;
A = 2;
N = [8,16,32,64];

% Initialize error vectors
rms_error = zeros(1,length(N));
max_error = zeros(1,length(N));

for ff=1:length(N)

    % FDS data
    if N(ff)==8;  M = dlmread('ns2d_8_devc.csv',',');  end
    if N(ff)==16; M = dlmread('ns2d_16_devc.csv',','); end
    if N(ff)==32; M = dlmread('ns2d_32_devc.csv',','); end
    if N(ff)==64; M = dlmread('ns2d_64_devc.csv',','); end

    t = M(:,1);
    U = M(:,2);
    P = M(:,3);
    nu = M(1,4);

    % grid
    dx(ff) = L/N(ff);
    dy(ff) = L/N(ff);
    x = dx(ff):dx(ff):L;
    y = dy(ff):dy(ff):L;
    xp = x-dx(ff)/2;
    yp = y-dy(ff)/2;

    % Exact solution
    ii = N(ff)/2;
    jj = N(ff)/2;
    u = 1 - A*cos(x(ii)-t).*sin(yp(jj)-t).*exp(-2*nu*t);
    p = -(A^2)/4*( cos(2*(xp(ii)-t)) + cos(2*(yp(jj)-t)) ).*exp(-4*nu*t);

    %figure(ff); plot(t,U,'bo'); hold on
    %figure(ff); plot(t,u,'b-')
    
    rms_error(ff) = sqrt( mean( (U-u).^2 ) );
    max_error(ff) = max(abs(U-u));

end

H(1)=loglog(dx,rms_error,'k*-'); hold on
H(2)=loglog(dx,dx,'k-');
H(3)=loglog(dx,dx.^2,'k--');
axis tight
xlabel('dx','fontsize',16)
ylabel('u-velocity rms error','fontsize',16)
legend(H(2:3),'O(dx)','O(dx^2)','Location','Northwest')

% write error file ------
M = [dx(1), dx(1), dx(1)^2, rms_error(1)];
for j = 2:length(dx)
    M = [ M; [dx(j), dx(j), dx(j)^2, rms_error(j)] ];
end
dlmwrite('error_inviscid.csv',M,'precision',6)
% -----------------------

clear all

L = 2*pi;
A = 2;
N = [8,16,32,64];

% Initialize error vectors
rms_error = zeros(1,length(N));
max_error = zeros(1,length(N));

for ff=1:length(N)

    % FDS data
    if N(ff)==8;  M = dlmread('ns2d_8_nupt1_devc.csv',',');  end
    if N(ff)==16; M = dlmread('ns2d_16_nupt1_devc.csv',','); end
    if N(ff)==32; M = dlmread('ns2d_32_nupt1_devc.csv',','); end
    if N(ff)==64; M = dlmread('ns2d_64_nupt1_devc.csv',','); end

    t = M(:,1);
    U = M(:,2);
    P = M(:,3);
    nu = M(1,4);

    % grid
    dx(ff) = L/N(ff);
    dy(ff) = L/N(ff);
    x = dx(ff):dx(ff):L;
    y = dy(ff):dy(ff):L;
    xp = x-dx(ff)/2;
    yp = y-dy(ff)/2;

    % Exact solution
    ii = N(ff)/2;
    jj = N(ff)/2;
    u = 1 - A*cos(x(ii)-t).*sin(yp(jj)-t).*exp(-2*nu*t);
    p = -(A^2)/4*( cos(2*(xp(ii)-t)) + cos(2*(yp(jj)-t)) ).*exp(-4*nu*t);

    figure
    G(1)=plot(t,U,'ko'); hold on
    G(2)=plot(t,u,'k-');
    xlabel('time, sec','fontsize',16)
    ylabel('u-velocity, m/s','fontsize',16)
    legend(G,'FDS','Analytical Solution','Location','Northeast')
    
    if (N(ff)==64)
        % analytical solution file ------
        M = [t(1), u(1)];
        for j = 2:length(t)
            M = [ M; [t(j), u(j)] ];
        end
        dlmwrite(['ns2d_nupt1_uvel_exact_soln.csv'],M,'precision',9)
        % -----------------------
    end

    rms_error(ff) = sqrt( mean( (U-u).^2 ) );
    max_error(ff) = max(abs(U-u));

end

figure
H(1)=loglog(dx,rms_error,'k*-'); hold on
H(2)=loglog(dx,dx,'k-');
H(3)=loglog(dx,dx.^2,'k--');
axis tight
xlabel('dx','fontsize',16)
ylabel('u-velocity rms error','fontsize',16)
legend(H(2:3),'O(dx)','O(dx^2)','Location','Northwest')

% write error file ------
M = [dx(1), dx(1), dx(1)^2, rms_error(1)];
for j = 2:length(dx)
    M = [ M; [dx(j), dx(j), dx(j)^2, rms_error(j)] ];
end
dlmwrite('error_viscous.csv',M,'precision',6)
% -----------------------

