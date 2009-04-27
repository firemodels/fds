% McDermott
% 2-19-09
% channel_flow.m

close all
clear all

M = dlmread('channel_devc.csv',',');
t = M(:,1);
ubar = M(:,2);
umax = M(:,3);
mu = max(M(:,4))
%rho = max(M(:,5))

plot(t,ubar,'linewidth',2)         % confirm steady state
xlabel('time (s)')
ylabel('velocity (m/s)')

rho = 1.198;               % fluid density
%mu = 1.81e-5;            % dynamic viscosity
nu = mu/rho;             % kinematic viscosity
H = 1;                   % channel height
dpdx = -100;               % mean pressure gradient (specified)
U = ubar(length(ubar))   % steady state mean velocity (planar averaged)

Re_H = U*H/nu                  % Reynolds number based on H
f_fds = 2*(-dpdx)*H/(rho*U^2)  % f from FDS
%f_exact = 24/Re_H              % exact value for Poiseuille flow

tau_w = H/4*abs(dpdx);
u_tau = sqrt(tau_w/rho);
z = 0.5*H./[8 16 32];
zp = z*u_tau/nu



