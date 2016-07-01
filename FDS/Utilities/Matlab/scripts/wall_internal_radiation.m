%Roger Wang
%6-29-11
%wall_internal_radiation.m

close all
clear all

infile  = '../../Verification/Radiation/wall_internal_radiation_devc.csv';
if ~exist(infile)
    display(['Error: File ',infile,' does not exist. Skipping case.'])
    return
end
M = csvread(infile,3,0);
   t = M(1,1);
   flux(1:5) = M(1,2:6)';
   
filename = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/wall_internal_radiation.tex';
fid = fopen(filename,'wt');

fprintf(fid,'%s\n','\begin{center}');
fprintf(fid,'%s\n', '\begin{tabular}{|c|c|c|} \hline');
fprintf(fid,'%s\n', '$\tau$      & $S(\tau)$   & FDS \\');
fprintf(fid,'%s\n', '            & (kW/m$^2$)  & (kW/m$^2$) \\ \hline\hline');
fprintf(fid,'%s %6.3f %s\n', '0.01        & 2.897       &',-flux(1),' \\');
fprintf(fid,'%s %5.2f %s\n', '0.1         & 24.94       &',-flux(2),' \\');
fprintf(fid,'%s %5.2f %s\n', '0.5         & 82.95       &',-flux(3),' \\');
fprintf(fid,'%s %5.1f %s\n', '1.0         & 116.3       &',-flux(4),' \\');
fprintf(fid,'%s %5.1f %s\n', '10.         & 149.0       &',-flux(5),' \\ \hline');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');

