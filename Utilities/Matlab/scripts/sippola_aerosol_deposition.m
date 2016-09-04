% Overholt
% 10-30-2012
% aerosol_depo.m
%
% Calculations for Aerosol Deposition FDS Validation cases
% located in (/Validation/Sippola_Aerosol_Deposition)

clear all

outdir = '../../../out/Sippola_Aerosol_Deposition/FDS_Output_Files/';

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
for i=1:length(filename)
    if ~exist([outdir, filename{i}])
        display(['Error: File ',[outdir, filename{i}],' does not exist. Skipping case.'])
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
            
% Particle Diameter (m)
particle_diameter = {1.0e-6, 2.8e-6, 5.2e-6, 9.1e-6 ...
                     16.e-6, 1.0e-6, 1.0e-6, 3.1e-6 ...
                     5.2e-6, 9.8e-6, 16.e-6, 1.0e-6 ...
                     3.1e-6, 5.4e-6, 8.7e-6, 15e-6};
                 
% Exp. deposition sampling area (m^2)
depo_area = 0.1 * 0.2;

% Simulation time (s)
time = 100;
                 
% Primary calculations
for i=1:length(filename)
    test = importdata([outdir, filename{i}],',',2);

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
    
    V_d_ceiling{i} = (J_sum_ceiling / (depo_area * time)) / (4 * avg_concentration); %#ok<*SAGROW>
    V_d_wall{i}    = (J_sum_wall    / (depo_area * time)) / (4 * avg_concentration);
    V_d_floor{i}   = (J_sum_floor   / (depo_area * time)) / (4 * avg_concentration);
    
    V_d_ceiling_nondim{i} = V_d_ceiling{i} / friction_velocity{i};
    V_d_wall_nondim{i}    = V_d_wall{i}    / friction_velocity{i};
    V_d_floor_nondim{i}   = V_d_floor{i}   / friction_velocity{i};
end

% Write out results file
fid = fopen([outdir,'Sippola_All_Tests.csv'],'wt','n');
fprintf(fid,'Test,Air Velocity 1 um (m/s),Friction Velocity 1 um (m/s),Ceiling Deposition Velocity 1 um (m/s),Wall Deposition Velocity 1 um (m/s),Floor Deposition Velocity 1 um (m/s),Dimensionless Ceiling Deposition Velocity 1 um (m/s),Dimensionless Wall Deposition Velocity 1 um (m/s),Dimensionless Floor Deposition Velocity 1 um (m/s),Test,Air Velocity 3 um (m/s),Friction Velocity 3 um (m/s),Ceiling Deposition Velocity 3 um (m/s),Wall Deposition Velocity 3 um (m/s),Floor Deposition Velocity 3 um (m/s),Dimensionless Ceiling Deposition Velocity 3 um (m/s),Dimensionless Wall Deposition Velocity 3 um (m/s),Dimensionless Floor Deposition Velocity 3 um (m/s),Test,Air Velocity 5 um (m/s),Friction Velocity 5 um (m/s),Ceiling Deposition Velocity 5 um (m/s),Wall Deposition Velocity 5 um (m/s),Floor Deposition Velocity 5 um (m/s),Dimensionless Ceiling Deposition Velocity 5 um (m/s),Dimensionless Wall Deposition Velocity 5 um (m/s),Dimensionless Floor Deposition Velocity 5 um (m/s),Test,Air Velocity 9 um (m/s),Friction Velocity 9 um (m/s),Ceiling Deposition Velocity 9 um (m/s),Wall Deposition Velocity 9 um (m/s),Floor Deposition Velocity 9 um (m/s),Dimensionless Ceiling Deposition Velocity 9 um (m/s),Dimensionless Wall Deposition Velocity 9 um (m/s),Dimensionless Floor Deposition Velocity 9 um (m/s),Test,Air Velocity 16 um (m/s),Friction Velocity 16 um (m/s),Ceiling Deposition Velocity 16 um (m/s),Wall Deposition Velocity 16 um (m/s),Floor Deposition Velocity 16 um (m/s),Dimensionless Ceiling Deposition Velocity 16 um (m/s),Dimensionless Wall Deposition Velocity 16 um (m/s),Dimensionless Floor Deposition Velocity 16 um (m/s)\n');
fprintf(fid,'1,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 2,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 3,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 4,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 5,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{1}, friction_velocity{1}, V_d_ceiling{1}, V_d_wall{1}, V_d_floor{1}, V_d_ceiling_nondim{1}, V_d_wall_nondim{1}, V_d_floor_nondim{1}, air_velocity{2},  friction_velocity{2},  V_d_ceiling{2},  V_d_wall{2},  V_d_floor{2},  V_d_ceiling_nondim{2},  V_d_wall_nondim{2},  V_d_floor_nondim{2},  air_velocity{3},  friction_velocity{3},  V_d_ceiling{3},  V_d_wall{3},  V_d_floor{3},  V_d_ceiling_nondim{3},  V_d_wall_nondim{3},  V_d_floor_nondim{3},  air_velocity{4},  friction_velocity{4},  V_d_ceiling{4},  V_d_wall{4},  V_d_floor{4},  V_d_ceiling_nondim{4},  V_d_wall_nondim{4},  V_d_floor_nondim{4},  air_velocity{5},  friction_velocity{5},  V_d_ceiling{5},  V_d_wall{5},  V_d_floor{5},  V_d_ceiling_nondim{5},  V_d_wall_nondim{5},  V_d_floor_nondim{5});
fprintf(fid,'6,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 8,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 9,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 10, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 11, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{6}, friction_velocity{6}, V_d_ceiling{6}, V_d_wall{6}, V_d_floor{6}, V_d_ceiling_nondim{6}, V_d_wall_nondim{6}, V_d_floor_nondim{6}, air_velocity{8},  friction_velocity{8},  V_d_ceiling{8},  V_d_wall{8},  V_d_floor{8},  V_d_ceiling_nondim{8},  V_d_wall_nondim{8},  V_d_floor_nondim{8},  air_velocity{9},  friction_velocity{9},  V_d_ceiling{9},  V_d_wall{9},  V_d_floor{9},  V_d_ceiling_nondim{9},  V_d_wall_nondim{9},  V_d_floor_nondim{9},  air_velocity{10}, friction_velocity{10}, V_d_ceiling{10}, V_d_wall{10}, V_d_floor{10}, V_d_ceiling_nondim{10}, V_d_wall_nondim{10}, V_d_floor_nondim{10}, air_velocity{11}, friction_velocity{11}, V_d_ceiling{11}, V_d_wall{11}, V_d_floor{11}, V_d_ceiling_nondim{11}, V_d_wall_nondim{11}, V_d_floor_nondim{11});
fprintf(fid,'7,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 13, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 14, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 15, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 16, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{7}, friction_velocity{7}, V_d_ceiling{7}, V_d_wall{7}, V_d_floor{7}, V_d_ceiling_nondim{7}, V_d_wall_nondim{7}, V_d_floor_nondim{7}, air_velocity{13}, friction_velocity{13}, V_d_ceiling{13}, V_d_wall{13}, V_d_floor{13}, V_d_ceiling_nondim{13}, V_d_wall_nondim{13}, V_d_floor_nondim{13}, air_velocity{14}, friction_velocity{14}, V_d_ceiling{14}, V_d_wall{14}, V_d_floor{14}, V_d_ceiling_nondim{14}, V_d_wall_nondim{14}, V_d_floor_nondim{14}, air_velocity{15}, friction_velocity{15}, V_d_ceiling{15}, V_d_wall{15}, V_d_floor{15}, V_d_ceiling_nondim{15}, V_d_wall_nondim{15}, V_d_floor_nondim{15}, air_velocity{16}, friction_velocity{16}, V_d_ceiling{16}, V_d_wall{16}, V_d_floor{16}, V_d_ceiling_nondim{16}, V_d_wall_nondim{16}, V_d_floor_nondim{16});
fprintf(fid,'12, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{12}, friction_velocity{12}, V_d_ceiling{12}, V_d_wall{12}, V_d_floor{12}, V_d_ceiling_nondim{12}, V_d_wall_nondim{12}, V_d_floor_nondim{12});
fclose(fid);