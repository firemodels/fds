% K McGrattan
% 8-26-2013
% convective_cooling_convergence.m

close all
clear all

infile{1} = '../../Verification/Heat_Transfer/convective_cooling_p1_devc.csv';
infile{2} = '../../Verification/Heat_Transfer/convective_cooling_p05_devc.csv';
infile{3} = '../../Verification/Heat_Transfer/convective_cooling_p025_devc.csv';
infile{4} = '../../Verification/Heat_Transfer/convective_cooling_p01_devc.csv';
infile{5} = '../../Verification/Heat_Transfer/convective_cooling_p005_devc.csv';
infile{6} = '../../Verification/Heat_Transfer/convective_cooling_p0025_devc.csv';
infile{7} = '../../Verification/Heat_Transfer/convective_cooling_p00125_devc.csv';

dx(1) = 0.1;
dx(2) = 0.05;
dx(3) = 0.025;
dx(4) = 0.01;
dx(5) = 0.005;
dx(6) = 0.0025;
dx(7) = 0.00125;

for k=1:7
     if ~exist(infile{k})
         display(['Error: File ',infile{k},' does not exist. Skipping case.'])
         return
     end
     M_10 = csvread(infile{k}, 2, 0);
     T(k) = M_10(length(M_10(:,1)),2);
end

%Write the error files
dir = '../../Verification/Heat_Transfer/';

filename11 = [dir,'convective_cooling_error.csv'];
fid11 = fopen(filename11,'wt');

fprintf(fid11,'%s\n', 'dx,dx^2,relative error');
   
T_exact = 295.3011488157;
for j = 1:7
    fprintf(fid11,'%8.3e, %8.3e, %9.5e\n', 0.01*dx(j),0.1*dx(j)^2,-(T(j)-T_exact)/T_exact);
end
fclose(fid11);



