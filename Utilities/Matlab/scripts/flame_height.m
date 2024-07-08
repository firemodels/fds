% McDermott
% 2-4-11
% flame_height.m
%
% integrates HRRPUL(z) from *_line.csv file to determine L_f/D (normalized flame height)

check_hrr  % confirm heat release rate

outdir = '../../../out/Heskestad_Flame_Height/';
expdir = '../../../exp/Heskestad_Flame_Height/';

filename = {'Qs=p1_RI=05_devc.csv',    'Qs=p1_RI=10_devc.csv',    'Qs=p1_RI=20_devc.csv';   ...
            'Qs=p2_RI=05_devc.csv',    'Qs=p2_RI=10_devc.csv',    'Qs=p2_RI=20_devc.csv';   ...
            'Qs=p5_RI=05_devc.csv',    'Qs=p5_RI=10_devc.csv',    'Qs=p5_RI=20_devc.csv';   ...
            'Qs=1_RI=05_devc.csv',     'Qs=1_RI=10_devc.csv',     'Qs=1_RI=20_devc.csv';    ...
            'Qs=2_RI=05_devc.csv',     'Qs=2_RI=10_devc.csv',     'Qs=2_RI=20_devc.csv';    ...
            'Qs=5_RI=05_devc.csv',     'Qs=5_RI=10_devc.csv',     'Qs=5_RI=20_devc.csv';    ...
            'Qs=10_RI=05_devc.csv',    'Qs=10_RI=10_devc.csv',    'Qs=10_RI=20_devc.csv';   ...
            'Qs=20_RI=05_devc.csv',    'Qs=20_RI=10_devc.csv',    'Qs=20_RI=20_devc.csv';   ...
            'Qs=50_RI=05_devc.csv',    'Qs=50_RI=10_devc.csv',    'Qs=50_RI=20_devc.csv';   ...
            'Qs=100_RI=05_devc.csv',   'Qs=100_RI=10_devc.csv',   'Qs=100_RI=20_devc.csv';  ...
            'Qs=200_RI=05_devc.csv',   'Qs=200_RI=10_devc.csv',   'Qs=200_RI=20_devc.csv';  ...
            'Qs=500_RI=05_devc.csv',   'Qs=500_RI=10_devc.csv',   'Qs=500_RI=20_devc.csv';  ...
            'Qs=1000_RI=05_devc.csv',  'Qs=1000_RI=10_devc.csv',  'Qs=1000_RI=20_devc.csv'; ...
            'Qs=2000_RI=05_devc.csv',  'Qs=2000_RI=10_devc.csv',  'Qs=2000_RI=20_devc.csv'; ...
            'Qs=5000_RI=05_devc.csv',  'Qs=5000_RI=10_devc.csv',  'Qs=5000_RI=20_devc.csv'; ...
            'Qs=10000_RI=05_devc.csv', 'Qs=10000_RI=10_devc.csv', 'Qs=10000_RI=20_devc.csv'};

rho_inf = 1.2;
cp = 1;
T_inf = 293;
g = 9.81;
D = 1.13;
f=0.99;
%Q* =  .1  .2  .5    1    2    5    10    20    50    100    200    500    1000    2000    5000    10000
Qdot=[151 303 756 1513 3025 7564 15127 30255 75636 151273 302545 756363 1512725 3025450 7563625 15127250];

for i=1:16 % hrr loop
    for j=1:3 % resolution loop
        M = csvread([outdir,filename{i,j}],2,0);
        L(j) = M(end,4);  % 99th percentile
        Qstar = Qdot(i)/(rho_inf*cp*T_inf*sqrt(g)*D^(5/2));
    end
    W(i,1:4) = [Qstar L(1)/D L(2)/D L(3)/D];
end

fclose('all');

% Write a file with FDS-predicted flame heights

header1 = {'Q*','L/D (RI=5)','L/D (RI=10)','L/D (RI=20)'};
filename1 = [outdir,'FDS_Flame_Height.csv'];
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s\n',header1{:});
for i=1:16
    fprintf(fid,'%f, %f, %f, %f\n',W(i,:));
end
fclose(fid);

% Generate FDS results for Tamanini cases

filename{1} = [outdir,'FDS_Tamanini_RI=05.csv'];
filename{2} = [outdir,'FDS_Tamanini_RI=10.csv'];
filename{3} = [outdir,'FDS_Tamanini_RI=20.csv'];

fds_line_file = {'Qs=1500_RI=05_line.csv', 'Qs=p6_RI=05_line.csv', 'Qs=p3_RI=05_line.csv'; ...
                 'Qs=1500_RI=10_line.csv', 'Qs=p6_RI=10_line.csv', 'Qs=p3_RI=10_line.csv'; ...
                 'Qs=1500_RI=20_line.csv', 'Qs=p6_RI=20_line.csv', 'Qs=p3_RI=20_line.csv'};

Qstar = [1500 .6 .3];
header = {'z/L jet','Q jet','z/L 62','Q 62','z/L 31','Q 31'};

for j=1:3 % resolution loop

    clear z hrr

    fid = fopen(filename{j},'wt');

    fprintf(fid,'%s, %s, %s, %s, %s, %s\n',header{:});
    A = [];

    for k=1:3 % hrr loop

        M = importdata([outdir,fds_line_file{j,k}],',',2);
        z = M.data(:,1);
        dz = z(2)-z(1);
        hrrpul = M.data(:,2);
        Qdot_line = sum(hrrpul)*dz;
        f = 0.99;

        % determine flame height
        for n=1:length(z)
            hrr(n) = sum(hrrpul(1:n))*dz/Qdot_line; % cummulative heat release
        end
        kk = find(hrr>f,1);
        L(k) = z(kk-1)+dz*(f-hrr(kk-1))/(hrr(kk)-hrr(kk-1));

        A=[A,(z+dz/2)/L(k),hrr'];

    end % hrr loop

    for i=1:length(z)
        fprintf(fid,'%f, %f, %f, %f, %f, %f\n',A(i,:));
    end
    fclose(fid);
    
end % resolution loop

