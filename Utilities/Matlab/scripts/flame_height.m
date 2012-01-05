% McDermott
% 2-4-11
% flame_height.m
%
% replaces old flame_height.f90
%
% integrates HRRPUL(z) from *_line.csv file to determine L_f/D (normalized
% flame height)

%close all
%clear all

addpath('../../Validation/Heskestad_Flame_Height/FDS_Output_Files/');

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
f=0.99;
%Q* =  .1  .2  .5    1    2    5    10    20    50    100    200    500    1000    2000    5000    10000
Qdot=[151 303 756 1513 3025 7564 15127 30255 75636 151273 302545 756363 1512725 3025450 7563625 15127250];

for i=1:16 % hrr loop
    for j=1:3 % resolution loop
        
        M = csvread(filename{i,j},2,0);
        z = M(:,1); dz = z(2)-z(1);
        hrrpul = M(:,2);
        Qdot_line = sum(hrrpul)*dz;
        Qstar = Qdot(i)/(rho_inf*cp*T_inf*sqrt(g)*D^(5/2));
        LfD = 3.7*Qstar^(2/5)-1.02; % Heskestad correlation Lf/D
        
        % determine flame height
        for n=1:length(z)
            hrr(n) = sum(hrrpul(1:n))*dz; % cummulative heat release
        end
        k = find(hrr>f*Qdot_line,1);
        if (k>1) 
            L(j) = z(k-1)+dz*(f*Qdot(i)-hrr(k-1))/(hrr(k)-hrr(k-1));
        else
            L(j) = dz*f*Qdot(i)/hrr(k);
		end
        
    end % resolution loop
    
    W(i,1:4) = [Qstar L(1)/D L(2)/D L(3)/D];
    H(i,1:2) = [Qstar LfD];
    
end % hrr loop

fclose('all');

header1 = {'Q*','L/D (RI=5)','L/D (RI=10)','L/D (RI=20)'};
filename1 = '../../Validation/Heskestad_Flame_Height/FDS_Output_Files/FDS_Flame_Height.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s\n',header1{:});
for i=1:16
    fprintf(fid,'%f, %f, %f, %f\n',W(i,:));
end
fclose(fid);
  
header2 = {'Q*','L/D'};
filename2 = '../../Validation/Heskestad_Flame_Height/Experimental_Data/Heskestad_Correlation.csv';
fid = fopen(filename2,'wt');
fprintf(fid,'%s, %s\n',header2{:});
for i=1:16
    fprintf(fid,'%f, %f\n',H(i,:));
end
fclose(fid);

% Generate FDS results for Tamanini cases

close all
clear all

filename{1} = '../../Validation/Heskestad_Flame_Height/FDS_Output_Files/FDS_Tamanini_RI=05.csv';
filename{2} = '../../Validation/Heskestad_Flame_Height/FDS_Output_Files/FDS_Tamanini_RI=10.csv';
filename{3} = '../../Validation/Heskestad_Flame_Height/FDS_Output_Files/FDS_Tamanini_RI=20.csv';

fds_line_file = {'Qs=2000_RI=05_line.csv', 'Qs=p5_RI=05_line.csv', 'Qs=p2_RI=05_line.csv'; ...
	             'Qs=2000_RI=10_line.csv', 'Qs=p5_RI=10_line.csv', 'Qs=p2_RI=10_line.csv'; ...
				 'Qs=2000_RI=20_line.csv', 'Qs=p5_RI=20_line.csv', 'Qs=p2_RI=20_line.csv'};
			 
Qstar = [2000 .5 .2];
header = {'z/L jet','Q jet','z/L 62','Q 62','z/L 31','Q 31'};

for j=1:3 % resolution loop
	
	clear z hrr
	
	fid = fopen(filename{j},'wt');
	
	fprintf(fid,'%s, %s, %s, %s, %s, %s\n',header{:});
	A = [];
	
	for k=1:3 % hrr loop
		
		M = importdata(fds_line_file{j,k},',',2);
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




















