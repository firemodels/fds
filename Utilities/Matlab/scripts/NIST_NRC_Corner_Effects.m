% McGrattan
% 12-6-2017
% NIST_NRC_Corner_Effects.m
%
% Reads the _devc.csv file and writes output in a form appropriate for dataplot.m

outdir = '../../../out/NIST_NRC_Corner_Effects/FDS_Output_Files/';

casename{1}  = 'corner_200_kW';
casename{2}  = 'corner_300_kW';
casename{3}  = 'corner_400_kW';
casename{4}  = 'wall_200_kW';
casename{5}  = 'wall_300_kW';
casename{6}  = 'wall_400_kW';
casename{7}  = 'cabinet_01';
casename{8}  = 'cabinet_02';
casename{9}  = 'cabinet_03';
casename{10} = 'cabinet_04';
casename{11} = 'cabinet_05';
casename{12} = 'cabinet_06';
casename{13} = 'cabinet_07';
casename{14} = 'cabinet_08';
casename{15} = 'cabinet_09';
casename{16} = 'cabinet_10';
casename{17} = 'cabinet_11';
casename{18} = 'cabinet_12';

for j=1:18

M = importdata([outdir,casename{j},'_devc.csv'],',',2);

H = cell(2,4);
H(1,:) = {'s' 'C' 'C' 'C'};
H(2,:) = {'Time' 'Lower' 'Middle' 'Upper'};

fid = fopen([outdir,casename{j},'_Plume.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s\n',H{1,:});
fprintf(fid,'%s,%s,%s,%s\n',H{2,:});
n_times = length(M.data(:,1));
for i=1:n_times
    low_smooth  = mean(M.data(max(1,i-5):min(n_times,i+5),60:88));
    mid_smooth  = mean(M.data(max(1,i-5):min(n_times,i+5),31:59));
    high_smooth = mean(M.data(max(1,i-5):min(n_times,i+5), 2:30));
    low  = max(low_smooth);
    mid  = max(mid_smooth);
    high = max(high_smooth);
    %low =  max(M.data(i,60:88));
    %mid =  max(M.data(i,31:59));
    %high = max(M.data(i, 2:30));
    fprintf(fid,'%4.0f,%5.1f,%5.1f,%5.1f\n',M.data(i,1),low,mid,high);
end
fclose(fid);
clear M

end

