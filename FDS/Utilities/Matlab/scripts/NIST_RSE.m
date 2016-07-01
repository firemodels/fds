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

addpath('../../Validation/NIST_RSE_1994/FDS_Output_Files')

%------------------------------------------------
% Read in FDS output files
%------------------------------------------------

HRR= [50 75 100 150 200 300 400 500 600];

%---50kW---%
RSE_50_struct=importdata('NIST_RSE_1994_50_RI=5_devc.csv'); %FDS 50 kW
RSE_data_RI_5(:,:,1)=RSE_50_struct.data(end,:); % data
RSE_50_struct=importdata('NIST_RSE_1994_50_RI=10_devc.csv'); %FDS 50 kW
RSE_data_RI_10(:,:,1)=RSE_50_struct.data(end,:); % data
RSE_50_struct=importdata('NIST_RSE_1994_50_RI=20_devc.csv'); %FDS 50 kW
RSE_data_RI_20(:,:,1)=RSE_50_struct.data(end,:); % data

%---75kW---%
RSE_75_struct=importdata('NIST_RSE_1994_75_RI=5_devc.csv'); %FDS 75 kW
RSE_data_RI_5(:,:,2)=RSE_75_struct.data(end,:); % data
RSE_75_struct=importdata('NIST_RSE_1994_75_RI=10_devc.csv'); %FDS 75 kW
RSE_data_RI_10(:,:,2)=RSE_75_struct.data(end,:); % data
RSE_75_struct=importdata('NIST_RSE_1994_75_RI=20_devc.csv'); %FDS 75 kW
RSE_data_RI_20(:,:,2)=RSE_75_struct.data(end,:); % data

%---100kW---%
RSE_100_struct=importdata('NIST_RSE_1994_100_RI=5_devc.csv'); %FDS 100 kW
RSE_data_RI_5(:,:,3)=RSE_100_struct.data(end,:); % data
RSE_100_struct=importdata('NIST_RSE_1994_100_RI=10_devc.csv'); %FDS 100 kW
RSE_data_RI_10(:,:,3)=RSE_100_struct.data(end,:); % data
RSE_100_struct=importdata('NIST_RSE_1994_100_RI=20_devc.csv'); %FDS 100 kW
RSE_data_RI_20(:,:,3)=RSE_100_struct.data(end,:); % data

%---150kW---%
RSE_150_struct=importdata('NIST_RSE_1994_150_RI=5_devc.csv'); %FDS 150 kW
RSE_data_RI_5(:,:,4)=RSE_150_struct.data(end,:); % data
RSE_150_struct=importdata('NIST_RSE_1994_150_RI=10_devc.csv'); %FDS 150 kW
RSE_data_RI_10(:,:,4)=RSE_150_struct.data(end,:); % data
RSE_150_struct=importdata('NIST_RSE_1994_150_RI=20_devc.csv'); %FDS 150 kW
RSE_data_RI_20(:,:,4)=RSE_150_struct.data(end,:); % data

%---200kW---%
RSE_200_struct=importdata('NIST_RSE_1994_200_RI=5_devc.csv'); %FDS 200 kW
RSE_data_RI_5(:,:,5)=RSE_200_struct.data(end,:); % data
RSE_200_struct=importdata('NIST_RSE_1994_200_RI=10_devc.csv'); %FDS 200 kW
RSE_data_RI_10(:,:,5)=RSE_200_struct.data(end,:); % data
RSE_200_struct=importdata('NIST_RSE_1994_200_RI=20_devc.csv'); %FDS 200 kW
RSE_data_RI_20(:,:,5)=RSE_200_struct.data(end,:); % data

%---300kW---%
RSE_300_struct=importdata('NIST_RSE_1994_300_RI=5_devc.csv'); %FDS 300 kW
RSE_data_RI_5(:,:,6)=RSE_300_struct.data(end,:); % data
RSE_300_struct=importdata('NIST_RSE_1994_300_RI=10_devc.csv'); %FDS 300 kW
RSE_data_RI_10(:,:,6)=RSE_300_struct.data(end,:); % data
RSE_300_struct=importdata('NIST_RSE_1994_300_RI=20_devc.csv'); %FDS 300 kW
RSE_data_RI_20(:,:,6)=RSE_300_struct.data(end,:); % data

%---400kW---%
RSE_400_struct=importdata('NIST_RSE_1994_400_RI=5_devc.csv'); %FDS 400 kW
RSE_data_RI_5(:,:,7)=RSE_400_struct.data(end,:); % data
RSE_400_struct=importdata('NIST_RSE_1994_400_RI=10_devc.csv'); %FDS 400 kW
RSE_data_RI_10(:,:,7)=RSE_400_struct.data(end,:); % data
RSE_400_struct=importdata('NIST_RSE_1994_400_RI=20_devc.csv'); %FDS 400 kW
RSE_data_RI_20(:,:,7)=RSE_400_struct.data(end,:); % data

%---500kW---%
RSE_500_struct=importdata('NIST_RSE_1994_500_RI=5_devc.csv'); %FDS 500 kW
RSE_data_RI_5(:,:,8)=RSE_500_struct.data(end,:); % data
RSE_500_struct=importdata('NIST_RSE_1994_500_RI=10_devc.csv'); %FDS 500 kW
RSE_data_RI_10(:,:,8)=RSE_500_struct.data(end,:); % data
RSE_500_struct=importdata('NIST_RSE_1994_500_RI=20_devc.csv'); %FDS 500 kW
RSE_data_RI_20(:,:,8)=RSE_500_struct.data(end,:); % data

%---600kW---%
RSE_600_struct=importdata('NIST_RSE_1994_600_RI=5_devc.csv'); %FDS 600 kW
RSE_data_RI_5(:,:,9)=RSE_600_struct.data(end,:); % data
RSE_600_struct=importdata('NIST_RSE_1994_600_RI=10_devc.csv'); %FDS 600 kW
RSE_data_RI_10(:,:,9)=RSE_600_struct.data(end,:); % data
RSE_600_struct=importdata('NIST_RSE_1994_600_RI=20_devc.csv'); %FDS 600 kW
RSE_data_RI_20(:,:,9)=RSE_600_struct.data(end,:); % data

%------------------------------------------------
% Building FDS data file with average species/temp
% Function of HRR
%------------------------------------------------

for i=1:9
    RSE_Results(1,i)=HRR(i);
    for j=2:16
        RSE_Results(j,i)=RSE_data_RI_5(:,j,i);
        RSE_Results((length(RSE_data_RI_5))+j-1,i)=RSE_data_RI_10(:,j,i);
        RSE_Results(2*(length(RSE_data_RI_5)-1)+j,i)=RSE_data_RI_20(:,j,i);
    end
end

%------------------------------------------------
% FDS Data CSV File
%------------------------------------------------

header1 = {'HRR','O2Rear_FDS_RI_5','CO2Rear_FDS_RI_5','CORear_FDS_RI_5','UHRear_FDS_RI_5','H2ORear_FDS_RI_5',...
    'O2Front_FDS_RI_5','CO2Front_FDS_RI_5','COFront_FDS_RI_5','UHFront_FDS_RI_5','H2OFront_FDS_RI_5',...
    'TRSampA_FDS_RI_5','TRSampBB_FDS_RI_5','TFSampA_FDS_RI_5','TFSampBB_FDS_RI_5','ITER_RI=5','O2Rear_FDS_RI_10',...
    'CO2Rear_FDS_RI_10','CORear_FDS_RI_10','UHRear_FDS_RI_10','H2ORear_FDS_RI_10','O2Front_FDS_RI_10',...
    'CO2Front_FDS_RI_10','COFront_FDS_RI_10','UHFront_FDS_RI_10','H2OFront_FDS_RI_10','TRSampA_FDS_RI_10',...
    'TRSampBB_FDS_RI_10','TFSampA_FDS_RI_10','TFSampBB_FDS_RI_10','ITER_RI=10','O2Rear_FDS_RI_20','CO2Rear_FDS_RI_20',...
    'CORear_FDS_RI_20','UHRear_FDS_RI_20','H2ORear_FDS_RI_20','O2Front_FDS_RI_20','CO2Front_FDS_RI_20',...
    'COFront_FDS_RI_20','UHFront_FDS_RI_20','H2OFront_FDS_RI_20','TRSampA_FDS_RI_20','TRSampBB_FDS_RI_20',...
    'TFSampA_FDS_RI_20','TFSampBB_FDS_RI_20','ITER_RI=20'};
filename1 = '../../Validation/NIST_RSE_1994/FDS_Output_Files/NIST_RSE_1994_FDS.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',header1{:});
for i=1:9
        fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n ',RSE_Results(:,i));
end
fclose(fid);
