% McGrattan
% 10-29-2019
% Crown_Fires.m
%
% Read the Crown_Fires _devc.csv files and determine the rate of spread based on the time history of front position.
% Write the results to a file that will be plotted via dataplot.m

close all
clear all

outdir = '../../../out/Crown_Fires/';
file_name = dir([outdir,'*_devc.csv']);
n_files = length(file_name);

for i=1:n_files;

   M = importdata([outdir file_name(i).name],',',2);
   indices = find(700<=M.data(:,2) & M.data(:,2)<=900 & M.data(:,1)>30 & M.data(:,1)<300);

   wind_speed(i) = mean(M.data(indices,3));
   p = polyfit(M.data(indices,1),M.data(indices,2),1); 
   slope(i) = p(1);

   clear indices M

end

fid = fopen([outdir 'ROS.csv'],'wt','n');
fprintf(fid,'%s\n','km/h,m/min');
fprintf(fid,'%s\n','U,ROS');
for i=1:n_files
   fprintf(fid,'%4.1f,%6.2f\n',3.6*wind_speed(i),60*slope(i));
end
fclose(fid);

