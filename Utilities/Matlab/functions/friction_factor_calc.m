% McDermott
% 5-14-2009
% friction_factor_calc.m

function [f_fds,Re_H]=friction_factor_calc(dpdx,mu,H,filename)

M = csvread(filename,2,0);

ubar = M(:,2);
rho = max(M(:,5));

U = ubar(length(ubar));  % steady state mean velocity (planar averaged)
Re_H = H*U*rho/mu;              % Reynolds number based on H
f_fds = 2*(-dpdx)*H/(rho*U^2);  % f from FDS


