% McDermott
% 8-29-2011
% cat_mccaffrey.m
%
% Concatenates columns from coarse and fine line files

FDS_Output_Files = '../../Validation/McCaffrey_Plume/FDS_Output_Files/';

M1 = importdata([FDS_Output_Files,'McCaffrey_14_kW_line.csv'],',',2);
M2 = importdata([FDS_Output_Files,'McCaffrey_14_kW_coarse_line.csv'],',',2);
M3 = importdata([FDS_Output_Files,'McCaffrey_14_kW_fine_line.csv'],',',2);

H = cell(2,7);
H(1,:) = {'m' 'C' 'm/s' 'C' 'm/s' 'C' 'm/s'};
H(2,:) = [M1.textdata(2,1:3),M2.textdata(2,2:3),M3.textdata(2,2:3)];

D = [M1.data(:,1:3),M2.data(:,2:3),M3.data(:,2:3)];

fid = fopen([FDS_Output_Files,'McCaffrey_14_kW_line_cat.csv'],'wt','n');

fprintf(fid,'%s, %s, %s, %s, %s, %s, %s\n',H{1,:});
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s\n',H{2,:});
for i=1:numel(M1.data(:,1))
    fprintf(fid,'%f, %f, %f, %f, %f, %f, %f\n',D(i,:));
end
fclose(fid);