% McGrattan
% 11-3-2014
% combine_csiro.m
%
% Sew together fire front positions from multiple mesh simulation of grassland fire spread

outdir = '../../../out/CSIRO_Grassland_Fires/FDS_Output_Files/';

label = {'Case_C064','Case_F19'};

H = cell(2,2);
H(1,:) = {'s' 'm'};
H(2,:) = {'Time' 'Front'};

fid2 = fopen([outdir,'csiro_speeds.csv'],'wt','n');

for j=1:2

M = importdata([outdir,label{j},'_devc.csv'],',',2);

fid = fopen([outdir,label{j},'_devc_combined.csv'],'wt','n');

fprintf(fid,'%s, %s\n',H{1,:});
fprintf(fid,'%s, %s\n',H{2,:});
for i=1:numel(M.data(:,1))
  [yy,ii]=max(M.data(i,8:13));
  fprintf(fid,'%f, %f\n',[M.data(i,1) M.data(i,ii+1)]);
  t(i) = M.data(i,1);
  x(i) = M.data(i,ii+1);
end
speed = (x(50)-x(30))/(t(50)-t(30));
fprintf(fid2,'%s, %f\n',label{j},speed);

fclose(fid);

clear fid M

end

fclose(fid2);
