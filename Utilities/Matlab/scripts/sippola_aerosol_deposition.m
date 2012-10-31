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

% Friction velocity (m/s) for 16 test cases
friction_velocity = {0.12, 0.12, 0.12, 0.13, 0.12, 0.28, 0.26, 0.26 ...
                     0.27, 0.28, 0.28, 0.45, 0.42, 0.44, 0.46, 0.45};

% Air velocity (m/s)
air_velocity = {2.2, 2.2, 2.1, 2.2, 2.2, 5.3, 5.2, 5.2 ...
                5.4, 5.3, 5.3, 9.0, 9.0, 8.8, 9.2, 9.1};
                 
% Exp. deposition sampling area (m^2)
depo_area = 0.1 * 0.2;

% Simulation time (s)
time = 100;
                 
% Primary calculations
for i=1:16
    test = importdata([FDS_Output_Files, filename{i}],',',2);

    upstream_concentration   = test.data(end,2);
    downstream_concentration = test.data(end,3);
    avg_concentration = (upstream_concentration + downstream_concentration) / 2;
    
    J1_ceiling = test.data(end,4);
    J2_ceiling = test.data(end,5);
    J3_ceiling = test.data(end,6);
    J4_ceiling = test.data(end,7);
    
    J1_wall = test.data(end,8);
    J2_wall = test.data(end,9);
    J3_wall = test.data(end,10);
    J4_wall = test.data(end,11);
    
    J1_floor = test.data(end,12);
    J2_floor = test.data(end,13);
    J3_floor = test.data(end,14);
    J4_floor = test.data(end,15);
    
    J_sum_ceiling = J1_ceiling + J2_ceiling + J3_ceiling + J4_ceiling;
    J_sum_wall    = J1_wall    + J2_wall    + J3_wall    + J4_wall;
    J_sum_floor   = J1_floor   + J2_floor   + J3_floor   + J4_floor;
    
    V_d_ceiling = (J_sum_ceiling / (depo_area * time)) / (4 * avg_concentration);
    V_d_wall    = (J_sum_wall    / (depo_area * time)) / (4 * avg_concentration);
    V_d_floor   = (J_sum_floor   / (depo_area * time)) / (4 * avg_concentration);
    
    V_d_ceiling_nondim = V_d_ceiling / friction_velocity{i};
    V_d_wall_nondim    = V_d_wall    / friction_velocity{i};
    V_d_floor_nondim   = V_d_floor   / friction_velocity{i};
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