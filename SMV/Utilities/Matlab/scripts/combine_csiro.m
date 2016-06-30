% McGrattan
% 11-3-2014
% combine_csiro.m
%
% Sew together fire front positions from multiple mesh simulation of grassland fire spread

FDS_Output_Files = '../../FDS/Validation/CSIRO_Grassland_Fires/FDS_Output_Files/';

M1 = importdata([FDS_Output_Files,'Case_C064_devc.csv'],',',2);
M2 = importdata([FDS_Output_Files,'Case_F19_devc.csv'],',',2);

H = cell(2,2);
H(1,:) = {'s' 'm'};
H(2,:) = {'Time' 'Front'};

fid1 = fopen([FDS_Output_Files,'Case_C064_devc_combined.csv'],'wt','n');
fid2 = fopen([FDS_Output_Files,'Case_F19_devc_combined.csv'],'wt','n');

fprintf(fid1,'%s, %s\n',H{1,:});
fprintf(fid1,'%s, %s\n',H{2,:});
for i=1:numel(M1.data(:,1))
  [yy,ii]=max(M1.data(i,8:13));
  fprintf(fid1,'%f, %f\n',[M1.data(i,1) M1.data(i,ii+1)]);
end

fprintf(fid2,'%s, %s\n',H{1,:});
fprintf(fid2,'%s, %s\n',H{2,:});
for i=1:numel(M2.data(:,1))
  [yy,ii]=max(M2.data(i,8:13));
  fprintf(fid2,'%f, %f\n',[M2.data(i,1) M2.data(i,ii+1)]);
end

fclose(fid1);
fclose(fid2);

