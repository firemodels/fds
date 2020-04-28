% Overholt
% 10-30-2012
% aerosol_depo.m
%
% Calculations for Aerosol Deposition FDS Validation cases
% located in (/Validation/Sippola_Aerosol_Deposition)

clear all

outdir = '../../../out/Sippola_Aerosol_Deposition/';

filename = {'Sippola_Test_01_devc.csv', 'Sippola_Test_02_devc.csv' ...
            'Sippola_Test_03_devc.csv', 'Sippola_Test_04_devc.csv' ...
            'Sippola_Test_05_devc.csv', 'Sippola_Test_06_devc.csv' ...
            'Sippola_Test_07_devc.csv', 'Sippola_Test_08_devc.csv' ...
            'Sippola_Test_09_devc.csv', 'Sippola_Test_10_devc.csv' ...
            'Sippola_Test_11_devc.csv', 'Sippola_Test_12_devc.csv' ...
            'Sippola_Test_13_devc.csv', 'Sippola_Test_14_devc.csv' ...
            'Sippola_Test_15_devc.csv', 'Sippola_Test_16_devc.csv' ...
            'Sippola_Test_17_devc.csv', 'Sippola_Test_18_devc.csv' ...
            'Sippola_Test_19_devc.csv', 'Sippola_Test_20_devc.csv' ...
            'Sippola_Test_21_devc.csv', 'Sippola_Test_22_devc.csv' ...
            'Sippola_Test_23_devc.csv', 'Sippola_Test_24_devc.csv' ...
            'Sippola_Test_25_devc.csv', 'Sippola_Test_26_devc.csv' ...
            'Sippola_Test_27_devc.csv', 'Sippola_Test_28_devc.csv' ...
            'Sippola_Test_29_devc.csv', 'Sippola_Test_30_devc.csv' ...
            'Sippola_Test_31_devc.csv'};

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
                     0.27, 0.28, 0.28, 0.45, 0.42, 0.44, 0.46, 0.45 ...
                     0.16, 0.16, 0.16, 0.16, 0.16, 0.37, 0.37, 0.37 ...
                     0.38, 0.38, 0.62, 0.62, 0.62, 0.64, 0.64};

% Air velocity (m/s)
air_velocity = {2.2, 2.2, 2.1, 2.2, 2.2, 5.3, 5.2, 5.2 ...
                5.4, 5.3, 5.3, 9.0, 9.0, 8.8, 9.2, 9.1 ...
                2.2, 2.2, 2.2, 2.2, 2.2, 5.3, 5.2, 5.2 ...
                5.3, 5.3, 8.9, 8.7, 8.9, 8.9, 8.9};
            
% Particle Diameter (m)
particle_diameter = {1.0e-6, 2.8e-6, 5.2e-6, 9.1e-6 ...
                     16.e-6, 1.0e-6, 1.0e-6, 3.1e-6 ...
                     5.2e-6, 9.8e-6, 16.e-6, 1.0e-6 ...
                     3.1e-6, 5.4e-6, 8.7e-6, 15e-6 ...
                     1.0e-6, 3.0e-6, 5.3e-6, 8.4e-6 ...
                     13.e-6, 1.0e-6, 2.9e-6, 4.9e-6 ...
                     8.2e-6, 13.e-6, 1.0e-6, 2.8e-6 ...
                     5.0e-6, 8.4e-6, 13.e-6};
                 
               
% Primary calculations
for i=1:length(filename)
    test = importdata([outdir, filename{i}],',',2);
    
    V_d_ceiling{i} = test.data(end,7);
    V_d_wall{i}    = test.data(end,8);
    V_d_floor{i}   = test.data(end,6);
    
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
fprintf(fid,'12,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 18,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 19,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 20,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 21,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{12}, friction_velocity{12}, V_d_ceiling{12}, V_d_wall{12}, V_d_floor{12}, V_d_ceiling_nondim{12}, V_d_wall_nondim{12}, V_d_floor_nondim{12}, air_velocity{18},  friction_velocity{18},  V_d_ceiling{18},  V_d_wall{18},  V_d_floor{18},  V_d_ceiling_nondim{18},  V_d_wall_nondim{18},  V_d_floor_nondim{18},  air_velocity{19},  friction_velocity{19},  V_d_ceiling{19},  V_d_wall{19},  V_d_floor{19},  V_d_ceiling_nondim{19},  V_d_wall_nondim{19},  V_d_floor_nondim{19},  air_velocity{20},  friction_velocity{20},  V_d_ceiling{20},  V_d_wall{20},  V_d_floor{20},  V_d_ceiling_nondim{20},  V_d_wall_nondim{20},  V_d_floor_nondim{20},  air_velocity{21},  friction_velocity{21},  V_d_ceiling{21},  V_d_wall{21},  V_d_floor{21},  V_d_ceiling_nondim{21},  V_d_wall_nondim{21},  V_d_floor_nondim{21});
fprintf(fid,'17,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 23,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 24,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 25, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 26, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{17}, friction_velocity{17}, V_d_ceiling{17}, V_d_wall{17}, V_d_floor{17}, V_d_ceiling_nondim{17}, V_d_wall_nondim{17}, V_d_floor_nondim{17}, air_velocity{23},  friction_velocity{23},  V_d_ceiling{23},  V_d_wall{23},  V_d_floor{23},  V_d_ceiling_nondim{23},  V_d_wall_nondim{23},  V_d_floor_nondim{23},  air_velocity{24},  friction_velocity{24},  V_d_ceiling{24},  V_d_wall{24},  V_d_floor{24},  V_d_ceiling_nondim{24},  V_d_wall_nondim{24},  V_d_floor_nondim{24},  air_velocity{25}, friction_velocity{25}, V_d_ceiling{25}, V_d_wall{25}, V_d_floor{25}, V_d_ceiling_nondim{25}, V_d_wall_nondim{25}, V_d_floor_nondim{25}, air_velocity{26}, friction_velocity{26}, V_d_ceiling{26}, V_d_wall{26}, V_d_floor{26}, V_d_ceiling_nondim{26}, V_d_wall_nondim{26}, V_d_floor_nondim{26});
fprintf(fid,'22,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 28, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 29, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 30, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 31, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{22}, friction_velocity{22}, V_d_ceiling{22}, V_d_wall{22}, V_d_floor{22}, V_d_ceiling_nondim{22}, V_d_wall_nondim{22}, V_d_floor_nondim{22}, air_velocity{28}, friction_velocity{28}, V_d_ceiling{28}, V_d_wall{28}, V_d_floor{28}, V_d_ceiling_nondim{28}, V_d_wall_nondim{28}, V_d_floor_nondim{28}, air_velocity{29}, friction_velocity{29}, V_d_ceiling{29}, V_d_wall{29}, V_d_floor{29}, V_d_ceiling_nondim{29}, V_d_wall_nondim{29}, V_d_floor_nondim{29}, air_velocity{30}, friction_velocity{30}, V_d_ceiling{30}, V_d_wall{30}, V_d_floor{30}, V_d_ceiling_nondim{30}, V_d_wall_nondim{30}, V_d_floor_nondim{30}, air_velocity{31}, friction_velocity{31}, V_d_ceiling{31}, V_d_wall{31}, V_d_floor{31}, V_d_ceiling_nondim{31}, V_d_wall_nondim{31}, V_d_floor_nondim{31});
fprintf(fid,'27,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 2,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 3,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 4,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 5,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{1}, friction_velocity{1}, V_d_ceiling{1}, V_d_wall{1}, V_d_floor{1}, V_d_ceiling_nondim{1}, V_d_wall_nondim{1}, V_d_floor_nondim{1}, air_velocity{2},  friction_velocity{2},  V_d_ceiling{2},  V_d_wall{2},  V_d_floor{2},  V_d_ceiling_nondim{2},  V_d_wall_nondim{2},  V_d_floor_nondim{2},  air_velocity{3},  friction_velocity{3},  V_d_ceiling{3},  V_d_wall{3},  V_d_floor{3},  V_d_ceiling_nondim{3},  V_d_wall_nondim{3},  V_d_floor_nondim{3},  air_velocity{4},  friction_velocity{4},  V_d_ceiling{4},  V_d_wall{4},  V_d_floor{4},  V_d_ceiling_nondim{4},  V_d_wall_nondim{4},  V_d_floor_nondim{4},  air_velocity{5},  friction_velocity{5},  V_d_ceiling{5},  V_d_wall{5},  V_d_floor{5},  V_d_ceiling_nondim{5},  V_d_wall_nondim{5},  V_d_floor_nondim{5});
fprintf(fid,'1,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 8,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 9,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 10, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 11, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{6}, friction_velocity{6}, V_d_ceiling{6}, V_d_wall{6}, V_d_floor{6}, V_d_ceiling_nondim{6}, V_d_wall_nondim{6}, V_d_floor_nondim{6}, air_velocity{8},  friction_velocity{8},  V_d_ceiling{8},  V_d_wall{8},  V_d_floor{8},  V_d_ceiling_nondim{8},  V_d_wall_nondim{8},  V_d_floor_nondim{8},  air_velocity{9},  friction_velocity{9},  V_d_ceiling{9},  V_d_wall{9},  V_d_floor{9},  V_d_ceiling_nondim{9},  V_d_wall_nondim{9},  V_d_floor_nondim{9},  air_velocity{10}, friction_velocity{10}, V_d_ceiling{10}, V_d_wall{10}, V_d_floor{10}, V_d_ceiling_nondim{10}, V_d_wall_nondim{10}, V_d_floor_nondim{10}, air_velocity{11}, friction_velocity{11}, V_d_ceiling{11}, V_d_wall{11}, V_d_floor{11}, V_d_ceiling_nondim{11}, V_d_wall_nondim{11}, V_d_floor_nondim{11});
fprintf(fid,'6,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 13, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 14, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 15, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 16, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{7}, friction_velocity{7}, V_d_ceiling{7}, V_d_wall{7}, V_d_floor{7}, V_d_ceiling_nondim{7}, V_d_wall_nondim{7}, V_d_floor_nondim{7}, air_velocity{13}, friction_velocity{13}, V_d_ceiling{13}, V_d_wall{13}, V_d_floor{13}, V_d_ceiling_nondim{13}, V_d_wall_nondim{13}, V_d_floor_nondim{13}, air_velocity{14}, friction_velocity{14}, V_d_ceiling{14}, V_d_wall{14}, V_d_floor{14}, V_d_ceiling_nondim{14}, V_d_wall_nondim{14}, V_d_floor_nondim{14}, air_velocity{15}, friction_velocity{15}, V_d_ceiling{15}, V_d_wall{15}, V_d_floor{15}, V_d_ceiling_nondim{15}, V_d_wall_nondim{15}, V_d_floor_nondim{15}, air_velocity{16}, friction_velocity{16}, V_d_ceiling{16}, V_d_wall{16}, V_d_floor{16}, V_d_ceiling_nondim{16}, V_d_wall_nondim{16}, V_d_floor_nondim{16});
fprintf(fid,'7,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 18,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 19,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 20,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 21,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{12}, friction_velocity{12}, V_d_ceiling{12}, V_d_wall{12}, V_d_floor{12}, V_d_ceiling_nondim{12}, V_d_wall_nondim{12}, V_d_floor_nondim{12}, air_velocity{18},  friction_velocity{18},  V_d_ceiling{18},  V_d_wall{18},  V_d_floor{18},  V_d_ceiling_nondim{18},  V_d_wall_nondim{18},  V_d_floor_nondim{18},  air_velocity{19},  friction_velocity{19},  V_d_ceiling{19},  V_d_wall{19},  V_d_floor{19},  V_d_ceiling_nondim{19},  V_d_wall_nondim{19},  V_d_floor_nondim{19},  air_velocity{20},  friction_velocity{20},  V_d_ceiling{20},  V_d_wall{20},  V_d_floor{20},  V_d_ceiling_nondim{20},  V_d_wall_nondim{20},  V_d_floor_nondim{20},  air_velocity{21},  friction_velocity{21},  V_d_ceiling{21},  V_d_wall{21},  V_d_floor{21},  V_d_ceiling_nondim{21},  V_d_wall_nondim{21},  V_d_floor_nondim{21});
fprintf(fid,'12,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 23,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 24,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 25, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 26, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{17}, friction_velocity{17}, V_d_ceiling{17}, V_d_wall{17}, V_d_floor{17}, V_d_ceiling_nondim{17}, V_d_wall_nondim{17}, V_d_floor_nondim{17}, air_velocity{23},  friction_velocity{23},  V_d_ceiling{23},  V_d_wall{23},  V_d_floor{23},  V_d_ceiling_nondim{23},  V_d_wall_nondim{23},  V_d_floor_nondim{23},  air_velocity{24},  friction_velocity{24},  V_d_ceiling{24},  V_d_wall{24},  V_d_floor{24},  V_d_ceiling_nondim{24},  V_d_wall_nondim{24},  V_d_floor_nondim{24},  air_velocity{25}, friction_velocity{25}, V_d_ceiling{25}, V_d_wall{25}, V_d_floor{25}, V_d_ceiling_nondim{25}, V_d_wall_nondim{25}, V_d_floor_nondim{25}, air_velocity{26}, friction_velocity{26}, V_d_ceiling{26}, V_d_wall{26}, V_d_floor{26}, V_d_ceiling_nondim{26}, V_d_wall_nondim{26}, V_d_floor_nondim{26});
fprintf(fid,'17,  %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 28, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 29, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 30, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, 31, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{22}, friction_velocity{22}, V_d_ceiling{22}, V_d_wall{22}, V_d_floor{22}, V_d_ceiling_nondim{22}, V_d_wall_nondim{22}, V_d_floor_nondim{22}, air_velocity{28}, friction_velocity{28}, V_d_ceiling{28}, V_d_wall{28}, V_d_floor{28}, V_d_ceiling_nondim{28}, V_d_wall_nondim{28}, V_d_floor_nondim{28}, air_velocity{29}, friction_velocity{29}, V_d_ceiling{29}, V_d_wall{29}, V_d_floor{29}, V_d_ceiling_nondim{29}, V_d_wall_nondim{29}, V_d_floor_nondim{29}, air_velocity{30}, friction_velocity{30}, V_d_ceiling{30}, V_d_wall{30}, V_d_floor{30}, V_d_ceiling_nondim{30}, V_d_wall_nondim{30}, V_d_floor_nondim{30}, air_velocity{31}, friction_velocity{31}, V_d_ceiling{31}, V_d_wall{31}, V_d_floor{31}, V_d_ceiling_nondim{31}, V_d_wall_nondim{31}, V_d_floor_nondim{31});
fprintf(fid,'22, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{22}, friction_velocity{22}, V_d_ceiling{22}, V_d_wall{22}, V_d_floor{22}, V_d_ceiling_nondim{22}, V_d_wall_nondim{22}, V_d_floor_nondim{22});
fprintf(fid,'27, %f, %f, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e, %3.1e\n', air_velocity{27}, friction_velocity{27}, V_d_ceiling{27}, V_d_wall{27}, V_d_floor{27}, V_d_ceiling_nondim{27}, V_d_wall_nondim{27}, V_d_floor_nondim{27});
fclose(fid);