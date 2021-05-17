% McDermott
% 5-17-2021
% favre_test.m
%
% compare TEMPORAL_STATISTIC='FAVRE AVERAGE' with brute force Favre average

close all
clear all

% Output files

datadir = '../../Verification/Species/';
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';

% gather Favre averages from line file

L = importdata([datadir,'favre_test_line.csv'],',',2);
YL1 = L.data(1,find(strcmp(L.colheaders,'YTILDE_O2')));
YL2 = L.data(2,find(strcmp(L.colheaders,'YTILDE_O2')));
YL3 = L.data(3,find(strcmp(L.colheaders,'YTILDE_O2')));

% brute force integration

M = importdata([datadir,'favre_test_devc.csv'],',',2);
t = M.data(find(M.data(:,1)>5),1);
RHO1 = M.data(find(M.data(:,1)>5),find(strcmp(M.colheaders,'"RHO_1"')));
RHO2 = M.data(find(M.data(:,1)>5),find(strcmp(M.colheaders,'"RHO_2"')));
RHO3 = M.data(find(M.data(:,1)>5),find(strcmp(M.colheaders,'"RHO_3"')));
Y1 = M.data(find(M.data(:,1)>5),find(strcmp(M.colheaders,'"YO2_1"')));
Y2 = M.data(find(M.data(:,1)>5),find(strcmp(M.colheaders,'"YO2_2"')));
Y3 = M.data(find(M.data(:,1)>5),find(strcmp(M.colheaders,'"YO2_3"')));

NUM1 = 0;
NUM2 = 0;
NUM3 = 0;
DENOM1 = 0;
DENOM2 = 0;
DENOM3 = 0;
for i=2:length(t)
    dt = t(i)-t(i-1);
    NUM1 = NUM1 + RHO1(i)*Y1(i)*dt;
    NUM2 = NUM2 + RHO2(i)*Y2(i)*dt;
    NUM3 = NUM3 + RHO3(i)*Y3(i)*dt;
    DENOM1 = DENOM1 + RHO1(i)*dt;
    DENOM2 = DENOM2 + RHO2(i)*dt;
    DENOM3 = DENOM3 + RHO3(i)*dt;
end
YTILDE1 = NUM1/DENOM1;
YTILDE2 = NUM2/DENOM2;
YTILDE3 = NUM3/DENOM3;

% compute error and report if necessary

e1 = abs(YL1-YTILDE1);
e2 = abs(YL2-YTILDE2);
e3 = abs(YL3-YTILDE3);

tol = 1.E-4;
if e1>tol; disp(['Matlab Warning: e1 = ',num2str(e1),' in Species/favre_test']); end
if e2>tol; disp(['Matlab Warning: e2 = ',num2str(e2),' in Species/favre_test']); end
if e3>tol; disp(['Matlab Warning: e3 = ',num2str(e3),' in Species/favre_test']); end



















