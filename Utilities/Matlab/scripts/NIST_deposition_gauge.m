% Overholt
% 10-30-2012
% NIST_deposition_gauge.m
%
% Thermophoretic Deposition FDS Validation cases
% located in (/Validation/NIST_Deposition_Gauge)

clear all

outdir = '../../../out/NIST_Deposition_Gauge/';

filename = {'NIST_SDG-2p5-100-p39a_devc.csv', 'NIST_SDG-2p5-100-p39b_devc.csv', 'NIST_SDG-2p5-100-p39c_devc.csv', 'NIST_SDG-2p5-100-p39d_devc.csv' ...
            'NIST_SDG-5p0-100-p39a_devc.csv', 'NIST_SDG-5p0-100-p39b_devc.csv', 'NIST_SDG-5p0-100-p39c_devc.csv', 'NIST_SDG-5p0-100-p39d_devc.csv' ...
            'NIST_SDG-10p0-100-p39a_devc.csv', 'NIST_SDG-10p0-100-p39b_devc.csv', 'NIST_SDG-10p0-100-p39c_devc.csv' ...
            'NIST_SDG-2p5-200-p39a_devc.csv', 'NIST_SDG-2p5-200-p39b_devc.csv', 'NIST_SDG-2p5-200-p39c_devc.csv', 'NIST_SDG-2p5-200-p39d_devc.csv' ...
            'NIST_SDG-5p0-200-p39a_devc.csv', 'NIST_SDG-5p0-200-p39b_devc.csv', 'NIST_SDG-5p0-200-p39c_devc.csv', 'NIST_SDG-5p0-200-p39d_devc.csv' ...
            'NIST_SDG-10p0-200-p39a_devc.csv', 'NIST_SDG-10p0-200-p39b_devc.csv', 'NIST_SDG-10p0-200-p39c_devc.csv' };

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

slpm = {2.5,2.5,2.5,2.5,5,5,5,5,10,10,10,2.5,2.5,2.5,2.5,5,5,5,5,10,10,10};

for i=1:22
    test = importdata([outdir, filename{i}],',',2);

    sensor1{i} = test.data(end,3);    
    sensor2{i} = test.data(end,4);
    sensor3{i} = test.data(end,5);
    sensor4{i} = test.data(end,6);
end

% Write out results file
fid = fopen([outdir,'SDG_All_Tests_100.csv'],'wt','n');
fprintf(fid,'SLPM,Sensor1,Sensor2,Sensor3,Sensor4\n');
for i=1:11
   fprintf(fid,'%5.3e,%5.3e,%5.3e,%5.3e,%5.3e\n',slpm{i},sensor1{i},sensor2{i},sensor3{i},sensor4{i});
end
fclose(fid);

fid = fopen([outdir,'SDG_All_Tests_200.csv'],'wt','n');
fprintf(fid,'SLPM,Sensor1,Sensor2,Sensor3,Sensor4\n');
for i=1:11
   fprintf(fid,'%5.3e,%5.3e,%5.3e,%5.3e,%5.3e\n',slpm{i},sensor1{i},sensor2{i},sensor3{i},sensor4{i});
end
fclose(fid);
