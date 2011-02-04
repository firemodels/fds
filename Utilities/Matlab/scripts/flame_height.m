% McDermott
% 2-4-11
% flame_height.m
%
% replaces old flame_height.f90
%
% integrates HRRPUL(z) from *_line.csv file to determine L_f/D (normalized
% flame height)

close all
clear all

addpath('../../Validation/Flame_Height/FDS_Output_Files/');

% list of line files
filename = {'Qs=p1_RI=05_line.csv',    'Qs=p1_RI=10_line.csv',    'Qs=p1_RI=20_line.csv';   ...
            'Qs=p2_RI=05_line.csv',    'Qs=p2_RI=10_line.csv',    'Qs=p2_RI=20_line.csv';   ...
            'Qs=p5_RI=05_line.csv',    'Qs=p5_RI=10_line.csv',    'Qs=p5_RI=20_line.csv';   ...
            'Qs=1_RI=05_line.csv',     'Qs=1_RI=10_line.csv',     'Qs=1_RI=20_line.csv';    ...
            'Qs=2_RI=05_line.csv',     'Qs=2_RI=10_line.csv',     'Qs=2_RI=20_line.csv';    ...
            'Qs=5_RI=05_line.csv',     'Qs=5_RI=10_line.csv',     'Qs=5_RI=20_line.csv';    ...
            'Qs=10_RI=05_line.csv',    'Qs=10_RI=10_line.csv',    'Qs=10_RI=20_line.csv';   ...
            'Qs=20_RI=05_line.csv',    'Qs=20_RI=10_line.csv',    'Qs=20_RI=20_line.csv';   ...
            'Qs=50_RI=05_line.csv',    'Qs=50_RI=10_line.csv',    'Qs=50_RI=20_line.csv';   ...
            'Qs=100_RI=05_line.csv',   'Qs=100_RI=10_line.csv',   'Qs=100_RI=20_line.csv';  ...
            'Qs=200_RI=05_line.csv',   'Qs=200_RI=10_line.csv',   'Qs=200_RI=20_line.csv';  ...
            'Qs=500_RI=05_line.csv',   'Qs=500_RI=10_line.csv',   'Qs=500_RI=20_line.csv';  ...
            'Qs=1000_RI=05_line.csv',  'Qs=1000_RI=10_line.csv',  'Qs=1000_RI=20_line.csv'; ...
            'Qs=2000_RI=05_line.csv',  'Qs=2000_RI=10_line.csv',  'Qs=2000_RI=20_line.csv'; ...
            'Qs=5000_RI=05_line.csv',  'Qs=5000_RI=10_line.csv',  'Qs=5000_RI=20_line.csv'; ...
            'Qs=10000_RI=05_line.csv', 'Qs=10000_RI=10_line.csv', 'Qs=10000_RI=20_line.csv'};

rho_inf = 1.2;
cp = 1;
T_inf = 293;
g = 9.81;
D = 1.13;

for i=1:16 % hrr loop
    for j=1:3 % resolution loop
        
        M = csvread(filename{i,j},2,0);
        z = M(:,1); dz = z(2)-z(1);
        hrrpul = M(:,2);
        Qdot = sum(hrrpul)*dz;
        Qstar = Qdot/(rho_inf*cp*T_inf*sqrt(g)*D^(5/2));
        LfD = 3.7*Qstar^(2/5)-1.02; % Heskestad correlation Lf/D
        
        % determine flame height
        for n=1:length(z)
            hrr(n) = sum(hrrpul(1:n))*dz; % cummulative heat release
        end
        k = find(hrr>.99*Qdot,1);
        L(j) = z(k-1)+dz*(.99*Qdot-hrr(k-1))/(hrr(k)-hrr(k-1));
        
    end % resolution loop
    
    W(i,1:4) = [Qstar L(1)/D L(2)/D L(3)/D];
    H(i,1:2) = [Qstar LfD];
    
end % hrr loop

header1 = {'Q*','L/D (RI=5)','L/D (RI=10)','L/D (RI=20)'};
filename1 = '../../Validation/Flame_Height/FDS_Output_Files/FDS_Flame_Height.csv';
fid = fopen(filename1,'w');
fprintf(fid,'%s, %s, %s, %s\n',header1{:});
for i=1:16
    fprintf(fid,'%f, %f, %f, %f\n',W(i,:));
end
fclose(fid);
  
header2 = {'Q*','L/D'};
filename2 = '../../Validation/Flame_Height/Experimental_Data/Heskestad_Correlation.csv';
fid = fopen(filename2,'w');
fprintf(fid,'%s, %s\n',header2{:});
for i=1:16
    fprintf(fid,'%f, %f\n',H(i,:));
end
fclose(fid);



