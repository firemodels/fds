% McDermott
% 1/23/2022
% tunnel_linear_cp.m

close all
clear all

% cp = a*T+b
a = 0.1584;
b = 953.5650;

T1=293.15;
u1=1;
W=28.8663; % gives rho1=1.2000
P=101325;
R=8314.5;
rho1=P*W/(R*T1);

q=2.5133e8; % W/m3
dx=0.001; % m

% find T2 from energy balance
% requires quadratic formula; note here that (a/2) is the usual "a" of the quadratic formula

c = -(T1^2 * (a/2) + b*T1) - q*dx/(rho1*u1);
T2 = ( -b + sqrt(b^2-4*(a/2)*c) )/(2*(a/2))

rho2 = P*W/(R*T2)
u2 = rho1*u1/rho2

u_0 = 0.5*(u1+u2)

dp2 = -(u2^2-u_0^2)/2 * rho2
dp1 = -(u_0^2-u1^2)/2 * 0.5*(rho1+rho2)
dp = dp1+dp2

% compute jump in energy across the fire

cp1 = a*T1+b
cp2 = a*T2+b

dq = (u2*rho2*cp2*T2 - u1*rho1*cp1*T1)/dx