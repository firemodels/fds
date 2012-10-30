% Overholt
% 10-30-2012
% aerosol_depo.m
%
% Calculations for Aerosol Deposition FDS Validation cases
% located in (/Validation/Sippola_Aerosol_Deposition)

clear all

FDS_Output_Files = '../../Validation/Sippola_Aerosol_Deposition/FDS_Output_Files/';

filename = {'Sippola_Test_01_devc.csv', 'Sippola_Test_02_devc.csv' ...
            'Sippola_Test_03_devc.csv', 'Sippola_Test_04_devc.csv' ...
            'Sippola_Test_05_devc.csv', 'Sippola_Test_06_devc.csv' ...
            'Sippola_Test_07_devc.csv', 'Sippola_Test_08_devc.csv' ...
            'Sippola_Test_09_devc.csv', 'Sippola_Test_10_devc.csv' ...
            'Sippola_Test_11_devc.csv', 'Sippola_Test_12_devc.csv' ...
            'Sippola_Test_13_devc.csv', 'Sippola_Test_14_devc.csv' ...
            'Sippola_Test_15_devc.csv', 'Sippola_Test_16_devc.csv'};

% Check for missing files
skip_case = 0;
for i=1:16
    if ~exist([FDS_Output_Files, filename{i}])
        display(['Error: File ',[FDS_Output_Files, filename{i}],' does not exist. Skipping case.'])
        skip_case = 1;
    end
end
if skip_case
    return
end

% Primary calculations
for i=1:16
    test = importdata([FDS_Output_Files, filename{i}],',',2);

    upstream_concentration = test.data(end,2) * 1e6
    downstream_concentration = test.data(end,3) * 1e6
    avg_concentration = (upstream_concentration + downstream_concentration) / 2
end

% Write out results file
% fid = fopen([FDS_Output_Files,'Sippola_All_Cases.csv'],'wt','n');
% 
% fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H3{1,:});
% fprintf(fid,'%s, %s, %s, %s, %s, %s\n',H3{2,:});
% for i=1:numel(S1(:,1))
%     fprintf(fid,'%f, %f, %f, %f, %f, %f\n',D3(i,:));
% end
% fclose(fid);