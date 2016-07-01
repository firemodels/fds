% Roger Wang
% 6-28-11
% ns2d.m
% Process the output of the ns2d verification cases

close all
clear all

infile{1,1} = '../../Verification/NS_Analytical_Solution/ns2d_8_devc.csv';
infile{2,1} = '../../Verification/NS_Analytical_Solution/ns2d_16_devc.csv';
infile{3,1} = '../../Verification/NS_Analytical_Solution/ns2d_32_devc.csv';
infile{4,1} = '../../Verification/NS_Analytical_Solution/ns2d_64_devc.csv';
infile{1,2} = '../../Verification/NS_Analytical_Solution/ns2d_8_nupt1_devc.csv';
infile{2,2} = '../../Verification/NS_Analytical_Solution/ns2d_16_nupt1_devc.csv';
infile{3,2} = '../../Verification/NS_Analytical_Solution/ns2d_32_nupt1_devc.csv';
infile{4,2} = '../../Verification/NS_Analytical_Solution/ns2d_64_nupt1_devc.csv';

outfile{1,1} = '../../Verification/NS_Analytical_Solution/ns2d_8_exact.csv';
outfile{2,1} = '../../Verification/NS_Analytical_Solution/ns2d_16_exact.csv';
outfile{3,1} = '../../Verification/NS_Analytical_Solution/ns2d_32_exact.csv';
outfile{4,1} = '../../Verification/NS_Analytical_Solution/ns2d_64_exact.csv';
outfile{1,2} = '../../Verification/NS_Analytical_Solution/ns2d_8_nupt1_exact.csv';
outfile{2,2} = '../../Verification/NS_Analytical_Solution/ns2d_16_nupt1_exact.csv';
outfile{3,2} = '../../Verification/NS_Analytical_Solution/ns2d_32_nupt1_exact.csv';
outfile{4,2} = '../../Verification/NS_Analytical_Solution/ns2d_64_nupt1_exact.csv';

nu(1) = 0.0;
nu(2) = 0.1;

pi = 4.0*atan(1.0);
x = pi;
y = pi;

dx(1) = 2.0*pi/8.0;
dx(2) = 2.0*pi/16.0;
dx(3) = 2.0*pi/32.0;
dx(4) = 2.0*pi/64.0;

for  k=1:2
    for j=1:4
        M_11 = fopen(outfile{j,k}, 'wt');
        fprintf(M_11, '%s\n', 'Time, u-vel');
        if ~exist(infile{j,k})
            display(['Error: File ',infile{j,k},' does not exist. Skipping case.'])
            return
        end
        M_10 = csvread(infile{j,k}, 2, 0);
        rms(j,k) = 0.0;
        cnt(j,k) = 0.0;
        for i = 1:length(M_10(:,1))
            t(i) = M_10(i,1);
            u(i,j,k)= M_10(i,2);
            u_exact = 1.0 - 2.0*cos(x-t(i))*sin(y-0.5*dx(j)-t(i))*exp(-2.0*nu(k)*t(i));
            fprintf(M_11,'%8.3f, %9.5f\n', t(i), u_exact);
            rms(j,k) = rms(j,k) + (u(i,j,k)-u_exact).^2;
            cnt(j,k) = cnt(j,k) + 1.0;
        end
        fclose(M_11);
        rms(j,k) = sqrt(rms(j,k)/cnt(j,k));
    end
end

       
%Write the error files
dir = '../../Verification/NS_Analytical_Solution/';

filename11 = [dir,'ns2d_error.csv'];
fid11 = fopen(filename11,'wt');
filename12 = [dir,'ns2d_nupt1_error.csv'];
fid12 = fopen(filename12,'wt');

fprintf(fid11,'%s\n', 'dx,dx^2,rms error');
fprintf(fid12,'%s\n', 'dx,dx^2,rms error');
   
for j = 1:4
    fprintf(fid11,'%8.3f, %8.3f, %8.4f\n', dx(j),dx(j)^2, rms(j,1));
    fprintf(fid12,'%8.3f, %8.3f, %8.4f\n', dx(j),dx(j)^2, rms(j,2));
end
fclose(fid11);
fclose(fid12);



