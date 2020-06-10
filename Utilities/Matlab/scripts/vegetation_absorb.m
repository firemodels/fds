% McGrattan
% 6-1-2018
% vegetation_absorb.m
%

close all
clear all

outdir = '../../Verification/WUI/';
infile = 'vegetation_absorb_devc.csv';

if ~exist([outdir, infile])
    display(['Error: File ',infile,' does not exist. Skipping case.'])
    return
end
            
M = importdata([outdir,infile],',',2);

fid = fopen([outdir,'vegetation_absorb_FDS.csv'],'wt','n');
fprintf(fid,'%s\n','mpuv,rad');
mpuv = [0.0 0.1 0.2 0.4 0.8 1.6];
for i=1:6
   fprintf(fid,'%5.1f,%5.2f\n',mpuv(i),M.data(end,i+1));
end
fclose(fid);

