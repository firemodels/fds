% Roger Wang
% 6-28-11
% radiation_box.m

close all
clear all

x( 1) = 2.50E-02;
x( 2) = 7.50E-02;
x( 3) = 1.25E-01;
x( 4) = 1.75E-01;
x( 5) = 2.25E-01;
x( 6) = 2.75E-01;
x( 7) = 3.25E-01;
x( 8) = 3.75E-01;
x( 9) = 4.25E-01;
x(10) = 4.75E-01;
x(11) = 5.25E-01;
x(12) = 5.75E-01;
x(13) = 6.25E-01;
x(14) = 6.75E-01;
x(15) = 7.25E-01;
x(16) = 7.75E-01;
x(17) = 8.25E-01;
x(18) = 8.75E-01;
x(19) = 9.25E-01;
x(20) = 9.75E-01;

dir = '../../Verification/Radiation/';

infile{1}  = [dir,'radiation_box__20___50_devc.csv'];
infile{2}  = [dir,'radiation_box__20__100_devc.csv'];
infile{3}  = [dir,'radiation_box__20__300_devc.csv'];
infile{4}  = [dir,'radiation_box__20_1000_devc.csv'];
infile{5}  = [dir,'radiation_box__20_2000_devc.csv'];
infile{6}  = [dir,'radiation_box_100___50_devc.csv'];
infile{7}  = [dir,'radiation_box_100__100_devc.csv'];
infile{8}  = [dir,'radiation_box_100__300_devc.csv'];
infile{9}  = [dir,'radiation_box_100_1000_devc.csv'];
infile{10} = [dir,'radiation_box_100_2000_devc.csv'];

label{ 1} = 'Flux_20_50';
label{ 2} = 'Flux_20_100';
label{ 3} = 'Flux_20_300';
label{ 4} = 'Flux_20_1000';
label{ 5} = 'Flux_20_2000';
label{ 6} = 'Flux_100_50';
label{ 7} = 'Flux_100_100';
label{ 8} = 'Flux_100_300';
label{ 9} = 'Flux_100_1000';
label{10} = 'Flux_100_2000';

for n=1:10
   if ~exist(infile{n})
       display(['Error: File ',infile{n},' does not exist. Skipping case.'])
       return
   end
   M = csvread(infile{n},3,0);
   t = M(1,1);
   flux(1:20,n) = M(1,2:21)';
end


filename = [dir,'radiation_box_devc.csv'];
fid = fopen(filename,'wt');

fprintf(fid,'%s','Position, ');
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',label{:});
for i=1:20
   fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',x(i),flux(i,1:10));
end
fclose(fid);













