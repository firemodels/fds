% McDermott
% 11-1-2016
% water_evap_1_const_gamma.m
%
% Calculation of expected values for water_evap_1_const_gamma.fds verification case.
% See fds/Verification/Sprinklers_and_Sprays/water_evap_1_const_gamma.csv
%
% Time,Rel. Hum,h_gas,h_water,dens,temp,pres,vapor
% 8,2.109675497,-27660.3572,1.9757398E+001,0.01,153.8927169,-7902.9592,0.01
% 10,2.109675497,-27660.3572,1.9757398E+001,0.01,153.8927169,-7902.9592,0.01

close all
clear all
%format long g

H_L = -1.9757398E+004; % J/kg (see _devc file)
R = 8314.5; % Pa * m3 / (kmol * K)
gam = 1.4;
P_1 = 101325; % Pa
V = 1; % m3
T_1 = 200 + 273.15; % K
T_w = 20  + 273.15; % K
M_w = 0.01; % kg (known from problem statement, see input file)

W_a = 28.84834; % kg/kmol (see .out file)
W_w = 18.01528; % kg/kmol (see .out file)

% ideal gas specific heats
cv_a = R/W_a * 1/(gam-1);
cv_w = R/W_w * 1/(gam-1);
cp_a = R/W_a * gam/(gam-1);
cp_w = R/W_w * gam/(gam-1);

% compute initial mass of air
RHO_1 = P_1 * W_a / (R * T_1);
M_a = RHO_1 * V;

% determine final mass fraction of water vapor in gas
% and compute new mixture molecular weight
Y_w = M_w/(M_a+M_w);
Y_a = 1-Y_w;
W = 1/(Y_w/W_w + Y_a/W_a);

% compute final density
RHO_2 = (M_w+M_a)/V;

% form energy balance (internal energy remains constant) determine final temperature
% and from ideal gas determine final pressure
T_2 = ( M_a*cv_a*T_1 + H_L ) / (M_a*cv_a + M_w*cv_w);
P_2 = RHO_2*R*T_2/W;

dP = P_2-P_1; % change in pressure, Pa

dH = (M_a*cp_a + M_w*cp_w)*T_2 - M_a*cp_a*T_1; % change in gas enthalpy, kJ

% output the results to the working directory

ddir = '../../Verification/Sprinklers_and_Sprays/';

fid = fopen([ddir,'water_evap_1_const_gamma.csv'],'w+');

fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s\n','Time','Rel. Hum','h_gas','h_water','dens','temp','pres','vapor');

fprintf(fid,'%i,%11.9f,%14.7e,%14.7e,%6.4f,%10.7f,%14.7e,%6.4f\n', ...
    8,2.109675497,dH/1000,-H_L/1000,RHO_2-RHO_1,T_2-273.15,dP,M_w);

fprintf(fid,'%i,%11.9f,%14.7e,%14.7e,%6.4f,%10.7f,%14.7e,%6.4f\n', ...
    10,2.109675497,dH/1000,-H_L/1000,RHO_2-RHO_1,T_2-273.15,dP,M_w);

fclose(fid);
























