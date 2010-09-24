% McDermott
% 9-23-10
% plot_mass.m

function [H]=plot_mass(mass_file,flux,area,row,col)

M = csvread(mass_file,row,0);
t = M(:,1);
m = M(:,col);

H(1)=plot(t,flux*area*t,'k-'); hold on
H(2)=plot(t,m,'b--');
