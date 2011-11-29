%------------------------------------------------
% C Weinschenk
% Process FDS simulation data of NIST RSE 1994 Experiments
%------------------------------------------------

clear all;
close all;

%------------------------------------------------
% Adding paths for data files
% and post processing scripts
%------------------------------------------------

addpath('/../../Validation/NIST_RSE_1994/FDS_Output_Files')

%------------------------------------------------
% Read in FDS output files
%------------------------------------------------

HRR= [50 75 100 150 200 300 400 500 600];

RSE_50_struct=importdata('NIST_RSE_1994_50_devc.csv'); %FDS 50 kW
RSE_data(:,:,1)=RSE_50_struct.data; % data

RSE_75_struct=importdata('NIST_RSE_1994_75_devc.csv'); %FDS 75 kW
RSE_data(:,:,2)=RSE_75_struct.data; % data

RSE_100_struct=importdata('NIST_RSE_1994_100_devc.csv'); %FDS 100 kW
RSE_data(:,:,3)=RSE_100_struct.data; % data

RSE_150_struct=importdata('NIST_RSE_1994_150_devc.csv'); %FDS 150 kW
RSE_data(:,:,4)=RSE_150_struct.data; % data

RSE_200_struct=importdata('NIST_RSE_1994_200_devc.csv'); %FDS 200 kW
RSE_data(:,:,5)=RSE_200_struct.data; % data

RSE_300_struct=importdata('NIST_RSE_1994_300_devc.csv'); %FDS 300 kW
RSE_data(:,:,6)=RSE_300_struct.data; % data

RSE_400_struct=importdata('NIST_RSE_1994_400_devc.csv'); %FDS 400 kW
RSE_data(:,:,7)=RSE_400_struct.data; % data

RSE_500_struct=importdata('NIST_RSE_1994_500_devc.csv'); %FDS 500 kW
RSE_data(:,:,8)=RSE_500_struct.data; % data

RSE_600_struct=importdata('NIST_RSE_1994_600_devc.csv'); %FDS 600 kW
RSE_data(:,:,9)=RSE_600_struct.data; % data

%------------------------------------------------
% Building FDS data file with average species/temp
% Function of HRR
%------------------------------------------------
start_row=12; % row to begin averaging (steady-state)

for i=1:9
    for j=2:15
        RSE_Results(1,j,i)=mean(RSE_data(start_row:length(RSE_data),j,i));
        RSE_Results(1,1,i)=HRR(i);
    end
end

%------------------------------------------------
% FDS Data CSV File
%------------------------------------------------

header1 = {'HRR','O2Rear_FDS','CO2Rear_FDS','CORear_FDS','UHRear_FDS','H2ORear_FDS','O2Front_FDS','CO2Front_FDS',...
    'COFront_FDS','UHFront_FDS','H2ORear_FDS','TRSampA_FDS','TRSampBB_FDS','TFSampA_FDS','TFSampBB_FDS'};
filename1 = '../../Validation/NIST_RSE_1994/FDS_Output_Files/NIST_RSE_1994_FDS.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',header1{:});
for i=1:9
        fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n ',RSE_Results(1,:,i));
end
fclose(fid);


