% McGrattan
% 8-10-09
% reaction_rate.m
%
% input:
% T = temperature (K)
% Y = mass fraction
%
% output:
% r = reaction rate

function r = reaction_rate(T,Y)

global dTdt
global R0
global E
global A
global residue

r(1) = -A(1).*Y(1)*exp(-E(1)./(R0.*T))./dTdt;
r(2) = -A(2).*Y(2)*exp(-E(2)./(R0.*T))./dTdt;
r(3) = -residue(2)*r(2);
r=r';


