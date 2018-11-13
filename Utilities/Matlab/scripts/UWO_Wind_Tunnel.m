% McGrattan
% 08-28-2018
% UWO_Test7_Case1.m
%
% Reads and transposes line data at various times in simulation

outdir = '../../../out/UWO_Wind_Tunnel/FDS_Output_Files/';
expdir = '../../../exp/UWO_Wind_Tunnel/';

label = {'UWO_Test7_Case1_180_5p00mm',...
         'UWO_Test7_Case1_180_2p50mm',...
         'UWO_Test7_Case1_180_1p25mm',...
         'UWO_Test7_Case1_270_5p00mm',...
         'UWO_Test7_Case1_270_2p50mm',...
         'UWO_Test7_Case1_270_1p25mm'};

H = cell(2,13);
H(1,:) = {'x' 'Cp_mean_1' 'Cp_rms_1' 'Cp_min_1' 'Cp_max_1' 'Cp_mean_2' 'Cp_rms_2' 'Cp_min_2' 'Cp_max_2' 'Cp_mean_3' 'Cp_rms_3' 'Cp_min_3' 'Cp_max_3'};
H(2,:) = {'1.0' 'NaN' 'NaN' 'NaN' 'NaN' 'NaN' 'NaN'  'NaN' 'NaN' 'NaN' 'NaN' 'NaN' 'NaN'};

for j=1:6

height = 0.04;
if j<4 ; length=0.14 ; else ; length=0.09 ; end

if j==1 ; n_lee=8  ; n_roof=28  ; end
if j==2 ; n_lee=16 ; n_roof=56  ; end
if j==3 ; n_lee=32 ; n_roof=112 ; end
if j==4 ; n_lee=8  ; n_roof=18  ; end
if j==5 ; n_lee=16 ; n_roof=36  ; end
if j==6 ; n_lee=32 ; n_roof=72  ; end

M = importdata([outdir,label{j},'_line.csv'],',',2);

fid = fopen([outdir,label{j},'.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',H{1,:});
for i=1:n_lee
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',M.data(i,1),M.data(i,2),M.data(i,8),M.data(i,14),M.data(i,20),M.data(i,3),M.data(i,9),M.data(i,15),M.data(i,21),M.data(i,4),M.data(i,10),M.data(i,16),M.data(i,22));
end
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',H{2,:});
for i=1:n_roof
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',height+M.data(i,26)+length/2,M.data(i,27),M.data(i,31),M.data(i,35),M.data(i,39),M.data(i,28),M.data(i,32),M.data(i,36),M.data(i,40),M.data(i,29),M.data(i,33),M.data(i,37),M.data(i,41));
end
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n',H{2,:});
for i=n_lee:-1:1
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',2*height+length-M.data(i,1),M.data(i,5),M.data(i,11),M.data(i,17),M.data(i,23),M.data(i,6),M.data(i,12),M.data(i,18),M.data(i,24),M.data(i,7),M.data(i,13),M.data(i,19),M.data(i,25));
end

fclose(fid);

fid = fopen([outdir,label{j},'s.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s,%s\n',H{1,1:5});
for i=1:n_roof
   fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',M.data(i,26)+length/2,M.data(i,30),M.data(i,34),M.data(i,38),M.data(i,42));
end

fclose(fid);
clear M fid

end

