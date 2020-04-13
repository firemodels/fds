% Overholt
% 10-30-2012
% aerosol_depo.m
%
% Calculations for Aerosol Deposition FDS Validation cases
% located in (/Validation/Sippola_Aerosol_Deposition)

clear all

outdir = '../../../out/Ranz_Marshall/';

filename = {'Ranz_Marshall_Table_1_1_devc.csv', 'Ranz_Marshall_Table_1_2_devc.csv' ...
            'Ranz_Marshall_Table_1_3_devc.csv', 'Ranz_Marshall_Table_1_4_devc.csv' ...
            'Ranz_Marshall_Table_1_5_devc.csv', 'Ranz_Marshall_Table_1_6_devc.csv' ...
            'Ranz_Marshall_Table_1_7_devc.csv', 'Ranz_Marshall_Table_1_8_devc.csv' ...
            'Ranz_Marshall_Table_1_9_devc.csv', 'Ranz_Marshall_Table_1_10_devc.csv' ...
            'Ranz_Marshall_Table_1_11_devc.csv', 'Ranz_Marshall_Table_1_12_devc.csv' ...
            'Ranz_Marshall_Table_1_13_devc.csv', 'Ranz_Marshall_Table_1_14_devc.csv' ...
            'Ranz_Marshall_Table_1_15_devc.csv', 'Ranz_Marshall_Table_1_16_devc.csv' ...
            'Ranz_Marshall_Table_1_17_devc.csv', 'Ranz_Marshall_Table_1_18_devc.csv' ...
            'Ranz_Marshall_Table_1_19_devc.csv'};

% Check for missing files
skip_case = 0;
for i=1:length(filename)
    if ~exist([outdir, filename{i}])
        display(['Error: File ',[outdir, filename{i}],' does not exist. Skipping case.'])
        skip_case = 1;
    end
end
if skip_case
    return
end

for i=1:length(filename)
    test = importdata([outdir, filename{i}],',',2); 
    d2dt{i} = test.data(end,6);
end

% Write out results file
fid = fopen([outdir,'Ranz_Marshall_Table_1.csv'],'wt','n');
fprintf(fid,'Test,dD2dt\n');
fprintf(fid,'1,  %3.1e\n', d2dt{1});
fprintf(fid,'2,  %3.1e\n', d2dt{2});
fprintf(fid,'3,  %3.1e\n', d2dt{3});
fprintf(fid,'4,  %3.1e\n', d2dt{4});
fprintf(fid,'5,  %3.1e\n', d2dt{5});
fprintf(fid,'6,  %3.1e\n', d2dt{6});
fprintf(fid,'7,  %3.1e\n', d2dt{7});
fprintf(fid,'8,  %3.1e\n', d2dt{8});
fprintf(fid,'9,  %3.1e\n', d2dt{9});
fprintf(fid,'10,  %3.1e\n', d2dt{10});
fprintf(fid,'11,  %3.1e\n', d2dt{11});
fprintf(fid,'12,  %3.1e\n', d2dt{12});
fprintf(fid,'13,  %3.1e\n', d2dt{13});
fprintf(fid,'14,  %3.1e\n', d2dt{14});
fprintf(fid,'15,  %3.1e\n', d2dt{15});
fprintf(fid,'16,  %3.1e\n', d2dt{16});
fprintf(fid,'17,  %3.1e\n', d2dt{17});
fprintf(fid,'18,  %3.1e\n', d2dt{18});
fprintf(fid,'19,  %3.1e\n', d2dt{19});
fclose(fid);

filename = {'Ranz_Marshall_Table_2_1_devc.csv', 'Ranz_Marshall_Table_2_2_devc.csv' ...
            'Ranz_Marshall_Table_2_3_devc.csv', 'Ranz_Marshall_Table_2_4_devc.csv' ...
            'Ranz_Marshall_Table_2_5_devc.csv', 'Ranz_Marshall_Table_2_6_devc.csv' ...
            'Ranz_Marshall_Table_2_7_devc.csv', 'Ranz_Marshall_Table_2_8_devc.csv' ...
            'Ranz_Marshall_Table_2_9_devc.csv'};

% Check for missing files
skip_case = 0;
for i=1:length(filename)
    if ~exist([outdir, filename{i}])
        display(['Error: File ',[outdir, filename{i}],' does not exist. Skipping case.'])
        skip_case = 1;
    end
end
if skip_case
    return
end

for i=1:length(filename)
    test = importdata([outdir, filename{i}],',',2);    
    d2dt{i} = test.data(end,6);
end

% Write out results file
fid = fopen([outdir,'Ranz_Marshall_Table_2.csv'],'wt','n');
fprintf(fid,'Test,dD2dt\n');
fprintf(fid,'1,  %3.1e\n', d2dt{1});
fprintf(fid,'2,  %3.1e\n', d2dt{2});
fprintf(fid,'3,  %3.1e\n', d2dt{3});
fprintf(fid,'4,  %3.1e\n', d2dt{4});
fprintf(fid,'5,  %3.1e\n', d2dt{5});
fprintf(fid,'6,  %3.1e\n', d2dt{6});
fprintf(fid,'7,  %3.1e\n', d2dt{7});
fprintf(fid,'8,  %3.1e\n', d2dt{8});
fprintf(fid,'9,  %3.1e\n', d2dt{9});
fclose(fid);

filename = {'Ranz_Marshall_Table_3_1_devc.csv', 'Ranz_Marshall_Table_3_2_devc.csv' ...
            'Ranz_Marshall_Table_3_3_devc.csv', 'Ranz_Marshall_Table_3_4_devc.csv' ...
            'Ranz_Marshall_Table_3_5_devc.csv', 'Ranz_Marshall_Table_3_6_devc.csv' ...
            'Ranz_Marshall_Table_3_7_devc.csv', 'Ranz_Marshall_Table_3_8_devc.csv' ...
            'Ranz_Marshall_Table_3_9_devc.csv'};

% Check for missing files
skip_case = 0;
for i=1:length(filename)
    if ~exist([outdir, filename{i}])
        display(['Error: File ',[outdir, filename{i}],' does not exist. Skipping case.'])
        skip_case = 1;
    end
end
if skip_case
    return
end

for i=1:length(filename)
    test = importdata([outdir, filename{i}],',',2);    
    d2dt{i} = test.data(end,6);
end

% Write out results file
fid = fopen([outdir,'Ranz_Marshall_Table_3.csv'],'wt','n');
fprintf(fid,'Test,dD2dt\n');
fprintf(fid,'1,  %3.1e\n', d2dt{1});
fprintf(fid,'2,  %3.1e\n', d2dt{2});
fprintf(fid,'3,  %3.1e\n', d2dt{3});
fprintf(fid,'4,  %3.1e\n', d2dt{4});
fprintf(fid,'5,  %3.1e\n', d2dt{5});
fprintf(fid,'6,  %3.1e\n', d2dt{6});
fprintf(fid,'7,  %3.1e\n', d2dt{7});
fprintf(fid,'8,  %3.1e\n', d2dt{8});
fprintf(fid,'9,  %3.1e\n', d2dt{9});
fclose(fid);

filename = {'Ranz_Marshall_Table_4_1_devc.csv', 'Ranz_Marshall_Table_4_2_devc.csv' ...
            'Ranz_Marshall_Table_4_3_devc.csv', 'Ranz_Marshall_Table_4_4_devc.csv' ...
            'Ranz_Marshall_Table_4_5_devc.csv', 'Ranz_Marshall_Table_4_6_devc.csv' ...
            'Ranz_Marshall_Table_4_7_devc.csv', 'Ranz_Marshall_Table_4_8_devc.csv' ...
            'Ranz_Marshall_Table_4_9_devc.csv', 'Ranz_Marshall_Table_4_10_devc.csv' ...
            'Ranz_Marshall_Table_4_11_devc.csv', 'Ranz_Marshall_Table_4_12_devc.csv' ...
            'Ranz_Marshall_Table_4_13_devc.csv'};

% Check for missing files
skip_case = 0;
for i=1:length(filename)
    if ~exist([outdir, filename{i}])
        display(['Error: File ',[outdir, filename{i}],' does not exist. Skipping case.'])
        skip_case = 1;
    end
end
if skip_case
    return
end

for i=1:length(filename)
    test = importdata([outdir, filename{i}],',',2);    
    d2dt{i} = test.data(end,6);
end

% Write out results file
fid = fopen([outdir,'Ranz_Marshall_Table_4.csv'],'wt','n');
fprintf(fid,'Test,dD2dt\n');
fprintf(fid,'1,  %3.1e\n', d2dt{1});
fprintf(fid,'2,  %3.1e\n', d2dt{2});
fprintf(fid,'3,  %3.1e\n', d2dt{3});
fprintf(fid,'4,  %3.1e\n', d2dt{4});
fprintf(fid,'5,  %3.1e\n', d2dt{5});
fprintf(fid,'6,  %3.1e\n', d2dt{6});
fprintf(fid,'7,  %3.1e\n', d2dt{7});
fprintf(fid,'8,  %3.1e\n', d2dt{8});
fprintf(fid,'9,  %3.1e\n', d2dt{9});
fprintf(fid,'10,  %3.1e\n', d2dt{10});
fprintf(fid,'11,  %3.1e\n', d2dt{11});
fprintf(fid,'12,  %3.1e\n', d2dt{12});
fprintf(fid,'13,  %3.1e\n', d2dt{13});
fclose(fid);