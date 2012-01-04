% McDermott
% 3-31-11
% fds_simple_chemistry.m
%
% This script is a distilled version of the SIMPLE_CHEMISTRY routine in FDS and may be
% used as a check on reaction coefficients for primitive and lumped species.

close all
clear all

% define the fuel, C_m H_n O_a N_b
m = 1;
n = 4;
a = 0;
b = 0;

% defind hydrogen atomic fraction in soot

X_H = 0;

% define CO and soot yields

y_CO = 0;
y_s  = 0;

% define the element matrix (number of atoms [rows] for each primitive species [columns])

% 'my fuel','oxygen','nitrogen','water vapor','carbon dioxide','carbon monoxide','soot'

i_fuel            = 1;
i_oxygen          = 2;
i_nitrogen        = 3;
i_water_vapor     = 4;
i_carbon_dioxide  = 5;
i_carbon_monoxide = 6;
i_soot            = 7;

n_species = 7;

E = [m 0 0 0 1 1 (1-X_H); ... % C
     n 0 0 2 0 0 X_H;     ... % H
     a 2 0 1 2 1 0;       ... % O
     b 0 2 0 0 0 0]       ... % N
     
W = E'*[12.0107 1.0794 15.9994 14.0067]' % primitive species molecular weights

% define the volume fractions of the background and fuel

v_0 = [0 .21 .79 .0 .0 0 0];
v_0 = v_0/sum(v_0) % normalize

v_1 = [1 0 0 0 0 0 0];
v_1 = v_1/sum(v_1) % normalize

% the reaction coefficients for the product primitive species temporarily stored in v_2

v_2 = zeros(1,n_species);

% compute what we know so far

v_2(i_carbon_monoxide) = W(i_fuel)/W(i_carbon_monoxide)*y_CO
v_2(i_soot)            = W(i_fuel)/W(i_soot)*y_s

% linear system right hand side

b = E*(v_1' - v_2')

% matrix

L = [E*v_0',E(:,i_carbon_dioxide),E(:,i_water_vapor),E(:,i_nitrogen)]

% solve the system

x = inv(L)*b

nu_0                  = x(1) % background stoichiometric coefficient
v_2(i_carbon_dioxide) = x(2);
v_2(i_water_vapor)    = x(3);
v_2(i_nitrogen)       = x(4);

nu_1 = -1       % fuel stoich coeff
nu_2 = sum(v_2) % prod stoich coeff

v_2 = v_2/nu_2; % normalize volume fractions

% check mass balance (should be 0)

nu_0*sum(v_0*W) + nu_1*sum(v_1*W) + nu_2*sum(v_2*W)

% check primitive reaction coefficients

NU = nu_0*v_0 + nu_1*v_1 + nu_2*v_2;

% display fuel properties

Z2Y = [v_0',v_1',v_2']




