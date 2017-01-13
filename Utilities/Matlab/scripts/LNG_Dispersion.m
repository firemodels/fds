% McGrattan
% 12-6-2016
% LNG_Disperion.m
%
% Reads and transposes line data at various times in simulation

outdir = '../../../out/LNG_Dispersion/FDS_Output_Files/';
expdir = '../../../exp/LNG_Dispersion/';

label = {'Burro3','Burro7','Burro8','Burro9'...
         'Coyote3','Coyote5','Coyote6'...
         'Falcon1','Falcon3','Falcon4'...
         'MaplinSands27','MaplinSands34','MaplinSands35'};

c_start = {[ 2 23 44 59],[ 2 23 44 56],[ 2 23 44 59],[23 44 59]...
           [ 2 29 44 50 65],[ 2 29 44 53 68],[ 2 29 44 53 68]...
           [ 2 30 58],[ 2 30 58],[ 2 30 58]...
           [ 2  5  8 14 17 23 26 29],[ 2  5],[ 2  5  8]};
c_end   = {[22 43 58 70],[22 43 55 64],[22 43 58 70],[43 58 70]...
           [28 43 49 64 67],[28 43 52 67 73],[28 43 52 67 73]...
           [29 57 85],[29 57 85],[29 57 85]...
           [ 4  7 13 16 22 25 28 31],[ 4  7],[ 4  7 10]};

H = cell(2,2);
H(1,:) = {'m' '%'};
H(2,:) = {'x' 'X_CH4'};

for j=1:13

M = importdata([outdir,label{j},'_devc.csv'],',',2);
E = importdata([expdir,label{j},'_exp.csv'], ',',2);

fid = fopen([outdir,label{j},'.csv'],'wt','n');
fprintf(fid,'%s, %s\n',H{1,:});
fprintf(fid,'%s, %s\n',H{2,:});
for i=1:numel(E.data(:,1))
   fprintf(fid,'%5.1f, %6.2f\n',E.data(i,1),max(0.01,100.*max(max(M.data(:,c_start{j}(i):c_end{j}(i))))));
end
fclose(fid);

clear M E fid

end

