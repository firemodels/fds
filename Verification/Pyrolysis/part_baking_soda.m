% McDermott
% 3-30-22
% part_baking_soda.m

close all
clear all

rho = 2200;
R = 8.3145;
A = 3.4e11;
E = 103000;
T = 420;
k = A*exp(-E/(R*T))
r_0 = 2.5e-6;

t = linspace(0,10,101);

% t = linspace(0,1,101); % used for T=500 K

% r = (r_0^3 * exp(-k*t)).^(1/3); % first-order model

r = r_0*(1-k*t); % spherical contraction model

d = r*2e6
% m = 4/3*pi*r.^3 * rho * 1e9
