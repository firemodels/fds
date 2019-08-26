% McGrattan
% 12-6-2016
% LNG_Disperion.m
%
% Reads and transposes line data at various times in simulation

outdir = '../../../out/LNG_Dispersion/';
expdir = '../../../exp/LNG_Dispersion/';

label = {'Burro3','Burro7','Burro8','Burro9'...
         'Coyote3','Coyote5','Coyote6'...
         'Falcon1','Falcon3','Falcon4'...
         'MaplinSands27','MaplinSands34','MaplinSands35'};

c_start = {[ 2 20 32 41],[ 2 20 41 50],[ 2 20 41 62],[ 2 20 32]...
           [ 2 29 44 50 65],[ 2 29 44 53 68],[ 2 29 44 53 68]...
           [ 2 25 49],[ 2 21 44],[ 2 21 53]...
           [ 2  5  8 12 13 19 22 25],[ 2  5],[ 2  5  8]};
c_end   = {[19 31 40 46],[19 40 49 52],[19 40 61 73],[19 31 40]...
           [28 43 49 64 67],[28 43 52 67 73],[28 43 52 67 73]...
           [24 48 65],[20 43 61],[20 52 67]...
           [ 4  7 11 12 18 21 24 27],[ 4  7],[ 4  7 10]};

t_start = {[ 30  30  60 170],[ 40  40  80 140],[ 60 150 400 650],[ 40  80 140]...
           [ 50  50  50  50  50],[ 40  40  40  40  40],[ 38  38  38  38  38]...
           [100 125 150],[100 125 175],[100 150 150]...
           [  0   0   0   0   0   0   0   0],[ 0  0],[  0   0   0]};
t_end   = {[130 130 160 270],[180 180 220 280],[140 230 480 730],[180 220 280]...
           [100 100 100 100 100],[130 130 130 130 130],[108 108 108 108 108]...
           [200 225 250],[250 275 325],[350 400 400]...
           [160 160 160 160 160 160 160 160],[95 95],[135 135 135]};

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
   time = M.data(:,1);
   idx = find(time>t_start{j}(i) & time<t_end{j}(i));
   fprintf(fid,'%5.1f, %6.2f\n',E.data(i,1),max(0.01,100.*max(max(M.data(idx,c_start{j}(i):c_end{j}(i))))));
end
fclose(fid);

clear M E fid

end

