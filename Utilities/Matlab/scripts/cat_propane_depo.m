% Overholt
% 6-13-2012
% cat_propane_depo.m
%
% Concatenates columns from Propane flame deposition FDS cases (/Verification/Species)

% Condensed phase aerosol

clear all

FDS_Output_Files = '../../Verification/Species/';

skip_case = 0;
if ~exist([FDS_Output_Files,'propane_flame_deposition_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'propane_flame_deposition_devc.csv'],' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([FDS_Output_Files,'propane_flame_deposition_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'propane_flame_deposition_none_devc.csv'],' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([FDS_Output_Files,'propane_flame_deposition_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'propane_flame_deposition_gravitational_devc.csv'],' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([FDS_Output_Files,'propane_flame_deposition_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'propane_flame_deposition_thermophoretic_devc.csv'],' does not exist. Skipping case.'])
    skip_case = 1;
end
if ~exist([FDS_Output_Files,'propane_flame_deposition_devc.csv'])
    display(['Error: File ',[FDS_Output_Files,'propane_flame_deposition_turbulent_devc.csv'],' does not exist. Skipping case.'])
    skip_case = 1;
end
if skip_case
    return
end

M1 = importdata([FDS_Output_Files,'propane_flame_deposition_devc.csv'],',',2);
M2 = importdata([FDS_Output_Files,'propane_flame_deposition_none_devc.csv'],',',2);
M3 = importdata([FDS_Output_Files,'propane_flame_deposition_gravitational_devc.csv'],',',2);
M4 = importdata([FDS_Output_Files,'propane_flame_deposition_thermophoretic_devc.csv'],',',2);
M5 = importdata([FDS_Output_Files,'propane_flame_deposition_turbulent_devc.csv'],',',2);

H1 = cell(2,6);
H1(1,:) = {'s' 'kg' 'kg' 'kg' 'kg' 'kg'};
H1(2,:) = {'Time' 'depo_all' 'depo_none' 'depo_gravitational' 'depo_thermophoretic' 'depo_turbulent'};

D1 = [M1.data(:,:), M2.data(:,2), M3.data(:,2), M4.data(:,2), M5.data(:,2)];

fid = fopen([FDS_Output_Files,'propane_flame_deposition_cat_wall.csv'],'wt','n');

fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H1{1,:});
fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H1{2,:});
for i=1:numel(M1.data(:,1))
    fprintf(fid,'%f, %f, %f, %f, %f, %f\n',D1(i,:));
end
fclose(fid);

% Gas phase aerosol

N1 = importdata([FDS_Output_Files,'propane_flame_deposition_mass.csv'],',',2);
N2 = importdata([FDS_Output_Files,'propane_flame_deposition_none_mass.csv'],',',2);
N3 = importdata([FDS_Output_Files,'propane_flame_deposition_gravitational_mass.csv'],',',2);
N4 = importdata([FDS_Output_Files,'propane_flame_deposition_thermophoretic_mass.csv'],',',2);
N5 = importdata([FDS_Output_Files,'propane_flame_deposition_turbulent_mass.csv'],',',2);

H2 = cell(2,6);
H2(1,:) = {'s' 'kg' 'kg' 'kg' 'kg' 'kg'};
H2(2,:) = {'Time' 'depo_all' 'depo_none' 'depo_gravitational' 'depo_thermophoretic' 'depo_turbulent'};

D2 = [N1.data(:,1), N1.data(:,8), N2.data(:,8), N3.data(:,8), N4.data(:,8), N5.data(:,8)];

fid = fopen([FDS_Output_Files,'propane_flame_deposition_cat_gas.csv'],'wt','n');

fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H2{1,:});
fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H2{2,:});
for i=1:numel(N1.data(:,1))
    fprintf(fid,'%f, %f, %f, %f, %f, %f\n',D2(i,:));
end
fclose(fid);

% Total aerosol (sum of gas plus wall)

S1 = interp1(M1.data(:,1), M1.data(:,2), N1.data(:,1)) + N1.data(:,8);
S2 = interp1(M2.data(:,1), M2.data(:,2), N2.data(:,1)) + N2.data(:,8);
S3 = interp1(M3.data(:,1), M3.data(:,2), N3.data(:,1)) + N3.data(:,8);
S4 = interp1(M4.data(:,1), M4.data(:,2), N4.data(:,1)) + N4.data(:,8);
S5 = interp1(M5.data(:,1), M5.data(:,2), N5.data(:,1)) + N5.data(:,8);

H3 = cell(2,6);
H3(1,:) = {'s' 'kg' 'kg' 'kg' 'kg' 'kg'};
H3(2,:) = {'Time' 'depo_all' 'depo_none' 'depo_gravitational' 'depo_thermophoretic' 'depo_turbulent'};

D3 = [N1.data(:,1), S1, S2, S3, S4, S5];

fid = fopen([FDS_Output_Files,'propane_flame_deposition_cat_total.csv'],'wt','n');

fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H3{1,:});
fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H3{2,:});
for i=1:numel(S1(:,1))
    fprintf(fid,'%f, %f, %f, %f, %f, %f\n',D3(i,:));
end
fclose(fid);
