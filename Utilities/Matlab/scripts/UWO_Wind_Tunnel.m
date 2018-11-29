% McGrattan
% 11-21-2018
% UWO_Wind Tunnel.m
%
% Reads and transposes line data for UWO Wind Tunnel cases

outdir = '../../../out/UWO_Wind_Tunnel/FDS_Output_Files/';
expdir = '../../../exp/UWO_Wind_Tunnel/';

H = cell(2,5);
H(1,:) = {'x' 'Cp_mean' 'Cp_rms' 'Cp_min' 'Cp_max'};
H(2,:) = {'1.0' 'NaN' 'NaN' 'NaN' 'NaN'};

input_file = fopen([outdir,'UWO_inputs.csv']);
C = textscan(input_file,'%s %s %d8 %d8 %d8 %f %f %d8 %d8 %d8 %d8  %d8 %d8 %d8 %f %f %d8 %d8 %d8 %d8  %d8 %d8 %d8 %f %f %d8 %d8 %d8 %d8  %d8 %d8 %d8 %f %f %d8 %d8 %d8 %d8','Delimiter',',','Headerlines',1);
fclose(input_file);

FDS_line_file = C{1};
output_file = C{2};

for k=1:length(C{1});

   M = importdata([outdir,C{1}{k}],',',2);

   n_sides = 4;
   if C{30}(k)==0; n_sides=3; end
   if C{21}(k)==0; n_sides=2; end
   if C{12}(k)==0; n_sides=1; end

   fid = fopen([outdir,output_file{k}],'wt','n');
   fprintf(fid,'%s,%s,%s,%s,%s\n',H{1,:});

   dist = 0.;

   for j=1:n_sides

      coord_col  = C{ 3+(j-1)*9}(k);
      points     = C{ 4+(j-1)*9}(k);
      order_index= C{ 5+(j-1)*9}(k);
      x1         = C{ 6+(j-1)*9}(k);
      x2         = C{ 7+(j-1)*9}(k);
      mean_col   = C{ 8+(j-1)*9}(k);
      rms_col    = C{ 9+(j-1)*9}(k);
      min_col    = C{10+(j-1)*9}(k);
      max_col    = C{11+(j-1)*9}(k);

      if order_index==1;
         for i=1:points;
            fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',dist+M.data(i,coord_col)-x1,M.data(i,mean_col),M.data(i,rms_col),M.data(i,min_col),M.data(i,max_col));
         end
      else
         for i=points:-1:1;
            fprintf(fid,'%6.3f,%6.3f,%6.3f,%6.3f,%6.3f\n',dist+x1-M.data(i,coord_col),M.data(i,mean_col),M.data(i,rms_col),M.data(i,min_col),M.data(i,max_col));
         end
      end

      if j<n_sides; fprintf(fid,'%s,%s,%s,%s,%s\n',H{2,:}); end
    
      dist = dist + abs(x2-x1);
   end

   fclose(fid);

   clear M
end

