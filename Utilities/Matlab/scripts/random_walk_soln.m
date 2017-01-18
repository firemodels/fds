% McDermott
% 9/2/2014
% random_walk_soln.m
%
% This script must be run prior to dataplot both to generate the expected results
% and to normalize the PDPA results from the FDS line file.

close all
clear all

datadir = '../../Verification/WUI/';

z0=5;
nz=80;
L=10;
dz=L/nz;
z=linspace(dz/2,L-dz/2,nz);
zp=z-z0;

RHO = 1.199; % kg/m3
D_1 = 0.1/RHO; % 0.0834 m2/s
D_2 = 0.01/RHO; % 0.00834 m2/s
t = 15; % s

f_1 = 1/(4*pi*D_1*sqrt(t))*exp(-zp.^2/(4*D_1*t));
f_1 = f_1/sum(f_1*dz); % normalize

f_2 = 1/(4*pi*D_2*sqrt(t))*exp(-zp.^2/(4*D_2*t));
f_2 = f_2/sum(f_2*dz); % normalize

% sum(f_1*dz)
% sum(f_2*dz)
% return

% print expected results to data directory

fid = fopen([datadir,'random_walk.csv'],'wt','n');
fprintf(fid,'%s, %s, %s\n','z','f_1','f_2');
for i=1:numel(z)
    fprintf(fid,'%f, %f, %f\n',z(i),f_1(i),f_2(i));
end
fclose(fid);

% pdpa_radius = 0.5; % m
% pdpa_volume = 4/3*pi*pdpa_radius^3; % m^3

% normalize and overwrite the FDS results for random_walk_1

filename = [datadir,'random_walk_1_line.csv'];

if ~exist(filename) % skip_case_if

    display(['Error: File ' filename ' does not exist. Skipping case.'])

else

    M = importdata([datadir,'random_walk_1_line.csv'],',',2);
    npv_z = M.data(:,find(strcmp(M.colheaders,'npv-z')));
    dz = (max(npv_z)-min(npv_z))/(length(npv_z)-1);
    npv = M.data(:,find(strcmp(M.colheaders,'npv')));
    nt = sum(npv);

    fid = fopen([datadir,'random_walk_1.csv'],'wt','n');
    fprintf(fid,'%s, %s\n','m',' ');
    fprintf(fid,'%s, %s\n','npv-z','npv');
    for i=1:numel(npv_z)
        fprintf(fid,'%f, %f\n',npv_z(i),npv(i)/nt/dz);
    end
    fclose(fid);

    % figure(1)
    % plot(npv_z,npv/nt/dz,'k-')

end

% normalize and overwrite the FDS results for random_walk_2

filename = [datadir,'random_walk_2_line.csv'];

if ~exist(filename) % skip_case_if

    display(['Error: File ' filename ' does not exist. Skipping case.'])

else

    M = importdata([datadir,'random_walk_2_line.csv'],',',2);
    npv_z = M.data(:,find(strcmp(M.colheaders,'npv-z')));
    dz = (max(npv_z)-min(npv_z))/(length(npv_z)-1);
    npv = M.data(:,find(strcmp(M.colheaders,'npv')));
    nt = sum(npv);

    fid = fopen([datadir,'random_walk_2.csv'],'wt','n');
    fprintf(fid,'%s, %s\n','m',' ');
    fprintf(fid,'%s, %s\n','npv-z','npv');
    for i=1:numel(npv_z)
        fprintf(fid,'%f, %f\n',npv_z(i),npv(i)/nt/dz);
    end
    fclose(fid);

    % figure(2)
    % plot(npv_z,npv/nt/dz,'k-')

end








