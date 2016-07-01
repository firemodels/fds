%Roger Wang
%7-20-11
%ashrae_7.m

close all
clear all

infile{1} = '../../Verification/HVAC/ashrae_7_exp.csv';
infile{2} = '../../Verification/HVAC/ashrae7_fixed_flow_devc.csv';
infile{3} = '../../Verification/HVAC/ashrae7_quadratic_devc.csv';
infile{4} = '../../Verification/HVAC/ashrae7_table_devc.csv';

label{1} = 'Experiment';
label{2} = 'Fixed Flow & ';
label{3} = 'Quadratic & ';
label{4} = 'Table & ';

duct{1} = '1';
duct{2} = '2';
duct{3} = '3';
duct{4} = '4';
duct{5} = '5';
duct{6} = '56';
duct{7} = '6';
duct{8} = '7';

M = csvread(infile{1},6,1);

for n = 2:4
    if ~exist(infile{n})
        display(['Error: File ',infile{n},' does not exist. Skipping case.'])
        return
    end
    m = csvread(infile{n},2,1);
    pressure(n,1:8) = m(length(m(:,1)),1:8);
end

filename = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ashrae_7.tex';
fid = fopen(filename, 'wt');

fprintf(fid,'%s\n','\begin{center}');
fprintf(fid,'%s\n','\begin{tabular}{|c|c|c|c|c|c|c|c|c|} \hline');
fprintf(fid, '%s', 'Duct Number & ');
fprintf(fid, '%s\n', '1 & 2 & 3 & 4 & 5 & 56 & 6 & 7 \\ \hline');
fprintf(fid, '%s\n', 'Experiment &');                                                                       
fprintf(fid, '%6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f %s\n', M,'\\');
for i = 2:4
    fprintf(fid, '%s', label{i});
    fprintf(fid, '%6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f & %6.3f %s\n', pressure(i,1:8),'\\');
end
fprintf(fid,'%s\n','\hline');
fprintf(fid,'%s\n','\end{tabular}');
fprintf(fid,'%s\n','\end{center}');

    
    