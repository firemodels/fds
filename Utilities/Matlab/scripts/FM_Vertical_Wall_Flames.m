% McGrattan
% 3-15-2017
% FM_Vertical_Wall_Flames.m
%
% Reads the _devc.csv file and writes output in a form appropriate for dataplot.m

outdir = '../../../out/FM_Vertical_Wall_Flames/';

nts = 32;

M = importdata([outdir,'propylene_devc.csv'],',',2);

% Heat flux data

H = cell(2,4);
H(1,:) = {'m' 'kW/m2' 'kW/m2' 'kW/m2'};
H(2,:) = {'z' 'HF-14' 'HF-18' 'HF-20'};

fid = fopen([outdir,'propylene_hf.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s\n',H{1,:});
fprintf(fid,'%s,%s,%s,%s\n',H{2,:});
for i=1:40
   z = 0.05*i-0.025;
   fprintf(fid,'%5.3f,%6.2f,%6.2f,%6.2f\n',z,M.data([8 10 11],1+i));
end
fclose(fid);

clear fid H

% Temperature data

H = cell(2,6);
H(1,:) = {'mm' 'K'    'K'    'K'    'K'    'K'};
H(2,:) = {'x'  'T-10' 'T-12' 'T-14' 'T-18' 'T-20'};

fid  = fopen([outdir,'propylene_tc.csv'],'wt','n');
fid2 = fopen([outdir,'propylene_tmp.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s,%s,%s\n',H{1,:});
fprintf(fid,'%s,%s,%s,%s,%s,%s\n',H{2,:});
fprintf(fid2,'%s,%s,%s,%s,%s,%s\n',H{1,:});
fprintf(fid2,'%s,%s,%s,%s,%s,%s\n',H{2,:});
for i=1:50
   x = 3*i-1.5;
   fprintf(fid, '%5.3f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n',x,M.data([5 6 7 9 10],41+i)+273.15);
   fprintf(fid2,'%5.3f,%6.2f,%6.2f,%6.2f,%6.2f,%6.2f\n',x,M.data([5 6 7 9 10],345+i)+273.15);
end
fclose(fid);
fclose(fid2);

clear fid fid2 H

% Soot data

H = cell(2,6);
H(1,:) = {'g/m2/s' 'mm' 'mm' 'mm' 'mm' 'mm'};
H(2,:) = {'mdot'   '365 mm'  '527 mm' '771 mm' '1022 mm' '1317 mm'};

threshold = 0.0025;
fid = fopen([outdir,'propylene_soot.csv'],'wt','n');
fprintf(fid,'%s,%s,%s,%s,%s,%s\n',H{1,:});
fprintf(fid,'%s,%s,%s,%s,%s,%s\n',H{2,:});
for i=1:nts
   for j=1:5
      for k=1:49
         if (M.data(i+1,92+50*(j-1)+k-1)>=threshold) && (M.data(i+1,92+50*(j-1)+k)<threshold) 
            index = k;
         end
      end
      delta(j) = index*3-1.5;
   end
   mdot = 2*i;
   fprintf(fid,'%5.3f,%5.1f,%5.1f,%5.1f,%5.1f,%5.1f\n',mdot,delta(1:5));
end
fclose(fid);

clear fid

