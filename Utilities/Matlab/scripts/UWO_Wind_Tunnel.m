% McGrattan
% 08-28-2018
% UWO_test7_case1.m
%
% Reads and transposes line data at various times in simulation

outdir = '../../../out/UWO_Wind_Tunnel/FDS_Output_Files/';
expdir = '../../../exp/UWO_Wind_Tunnel/';

label = {'UWO_test7_case1_180_32',...
         'UWO_test7_case1_180_64',...
         'UWO_test7_case1_180_128',...
         'UWO_test7_case1_270_32',...
         'UWO_test7_case1_270_64',...
         'UWO_test7_case1_270_128'};

H = cell(2,7);
H(1,:) = {'x' 'Cp_mean_1' 'Cp_rms_1' 'Cp_mean_2' 'Cp_rms_2' 'Cp_mean_3' 'Cp_rms_3'};
H(2,:) = {'1.0' 'NaN' 'NaN' 'NaN' 'NaN' 'NaN' 'NaN'};

for j=1:6

height = 0.04;
if j<4 ; length=0.14 ; else ; length=0.09 ; end

M = importdata([outdir,label{j},'_line.csv'],',',2);

fid = fopen([outdir,label{j},'.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',H{1,:});
for i=1:10
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',M.data(i,1),M.data(i,2),M.data(i,4),M.data(i,14),M.data(i,16),M.data(i,26),M.data(i,28));
end
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',H{2,:});
for i=1:20
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',height+M.data(i,5)+length/2,M.data(i,6),M.data(i,8),M.data(i,18),M.data(i,20),M.data(i,30),M.data(i,32));
end
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s\n',H{2,:});
for i=10:-1:1
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',2*height+length-M.data(i,9),M.data(i,10),M.data(i,12),M.data(i,22),M.data(i,24),M.data(i,34),M.data(i,36));
end

fclose(fid);

fid = fopen([outdir,label{j},'s.csv'],'wt','n');
fprintf(fid,'%s,%s,%s\n',H{1,1:3});
for i=1:20
   fprintf(fid,'%6.3f,%6.3f,%6.3f\n',M.data(i,37)+length/2,M.data(i,38),M.data(i,40));
end

fclose(fid);
clear M fid

end

