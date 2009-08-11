% McGrattan
% 8-10-09
% reaction_rate.m
%
% input:
% T = Temperature (K)
% Y = mass fraction
%
% output:
% r = Reaction rate
%

function r = reaction_rate(T,Y)

dTdt = 5./60.;
R0 = 8314.3;
T_0 = 300.+273.;
%r_0 = 0.001;
%A = 1.e13;
%E = -R0*T_0*log(r_0./A);
E=400e6;
A=(E./(R0.*T_0.*T_0))*dTdt*exp(E./(R0.*T_0));

r = -A.*Y*exp(-E./(R0.*T))./dTdt;

