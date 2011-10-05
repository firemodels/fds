%------------------------------------------------
% C Weinschenk
% FDS Simulation data of NIST RSE 1994 Experiments
%------------------------------------------------

% clear all;
% close all;

%------------------------------------------------
% Adding paths for data files
% and post processing scripts
%------------------------------------------------

addpath('../../Validation/NIST_RSE_1994/Test/')

%------------------------------------------------
% Read in FDS output files
% 9 FDS HRRs
%------------------------------------------------

%-------------
% 50 kW
%-------------
RSE_50_struct=importdata('NIST_RSE_1994_50_devc.csv'); %FDS data file
RSE_50_header=RSE_50_struct.textdata; % headers
RSE_50data=RSE_50_struct.data; % data

time50=RSE_50data(:,1); % time

O2_rear_fds_50=RSE_50data(:,2); % oxygen volume percent rear
CO2_rear_fds_50=RSE_50data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_50=RSE_50data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_50=RSE_50data(:,5); % fuel volume percent rear
H2O_rear_fds_50=RSE_50data(:,6); % water vapor volume percent rear

O2_front_fds_50=RSE_50data(:,7); % oxygen volume percent front
CO2_front_fds_50=RSE_50data(:,8); % carbon dioxide volume percent front
CO_front_fds_50=RSE_50data(:,9); % carbon monoxide volume percent front
Methane_front_fds_50=RSE_50data(:,10); % fuel volume percent front
H2O_front_fds_50=RSE_50data(:,11); % water vapor volume percent front

temp_rearA_fds_50=RSE_50data(:,12); % temperature rear sample A
temp_rearB_fds_50=RSE_50data(:,13); % temperature rear sample B

temp_frontA_fds_50=RSE_50data(:,14); % temperature front sample A
temp_frontB_fds_50=RSE_50data(:,15); % temperature front sample 
% -------------
% 75 kW
% -------------
RSE_75_struct=importdata('NIST_RSE_1994_75_devc.csv'); %FDS data file
RSE_75_header=RSE_75_struct.textdata; % headers
RSE_75data=RSE_75_struct.data; % data

time75=RSE_75data(:,1); % time

O2_rear_fds_75=RSE_75data(:,2); % oxygen volume percent rear
CO2_rear_fds_75=RSE_75data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_75=RSE_75data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_75=RSE_75data(:,5); % fuel volume percent rear
H2O_rear_fds_75=RSE_75data(:,6); % water vapor volume percent rear

O2_front_fds_75=RSE_75data(:,7); % oxygen volume percent front
CO2_front_fds_75=RSE_75data(:,8); % carbon dioxide volume percent front
CO_front_fds_75=RSE_75data(:,9); % carbon monoxide volume percent front
Methane_front_fds_75=RSE_75data(:,10); % fuel volume percent front
H2O_front_fds_75=RSE_75data(:,11); % water vapor volume percent front

temp_rearA_fds_75=RSE_75data(:,12); % temperature rear sample A
temp_rearB_fds_75=RSE_75data(:,13); % temperature rear sample B

temp_frontA_fds_75=RSE_75data(:,14); % temperature front sample A
temp_frontB_fds_75=RSE_75data(:,15); % temperature front sample B
% -------------
% 100 kW
% -------------
RSE_100_struct=importdata('NIST_RSE_1994_100_devc.csv'); %FDS data file
RSE_100_header=RSE_100_struct.textdata; % headers
RSE_100data=RSE_100_struct.data; % data

time100=RSE_100data(:,1); % time

O2_rear_fds_100=RSE_100data(:,2); % oxygen volume percent rear
CO2_rear_fds_100=RSE_100data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_100=RSE_100data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_100=RSE_100data(:,5); % fuel volume percent rear
H2O_rear_fds_100=RSE_100data(:,6); % water vapor volume percent rear

O2_front_fds_100=RSE_100data(:,7); % oxygen volume percent front
CO2_front_fds_100=RSE_100data(:,8); % carbon dioxide volume percent front
CO_front_fds_100=RSE_100data(:,9); % carbon monoxide volume percent front
Methane_front_fds_100=RSE_100data(:,10); % fuel volume percent front
H2O_front_fds_100=RSE_100data(:,11); % water vapor volume percent front

temp_rearA_fds_100=RSE_100data(:,12); % temperature rear sample A
temp_rearB_fds_100=RSE_100data(:,13); % temperature rear sample B

temp_frontA_fds_100=RSE_100data(:,14); % temperature front sample A
temp_frontB_fds_100=RSE_100data(:,15); % temperature front sample B
% -------------
% 150 kW
% -------------
RSE_150_struct=importdata('NIST_RSE_1994_150_devc.csv'); %FDS data file
RSE_150_header=RSE_150_struct.textdata; % headers
RSE_150data=RSE_150_struct.data; % data

time150=RSE_150data(:,1); % time

O2_rear_fds_150=RSE_150data(:,2); % oxygen volume percent rear
CO2_rear_fds_150=RSE_150data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_150=RSE_150data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_150=RSE_150data(:,5); % fuel volume percent rear
H2O_rear_fds_150=RSE_150data(:,6); % water vapor volume percent rear

O2_front_fds_150=RSE_150data(:,7); % oxygen volume percent front
CO2_front_fds_150=RSE_150data(:,8); % carbon dioxide volume percent front
CO_front_fds_150=RSE_150data(:,9); % carbon monoxide volume percent front
Methane_front_fds_150=RSE_150data(:,10); % fuel volume percent front
H2O_front_fds_150=RSE_150data(:,11); % water vapor volume percent front

temp_rearA_fds_150=RSE_150data(:,12); % temperature rear sample A
temp_rearB_fds_150=RSE_150data(:,13); % temperature rear sample B

temp_frontA_fds_150=RSE_150data(:,14); % temperature front sample A
temp_frontB_fds_150=RSE_150data(:,15); % temperature front sample B
% -------------
% 200 kW
% -------------
RSE_200_struct=importdata('NIST_RSE_1994_200_devc.csv'); %FDS data file
RSE_200_header=RSE_200_struct.textdata; % headers
RSE_200data=RSE_200_struct.data; % data

time200=RSE_200data(:,1); % time

O2_rear_fds_200=RSE_200data(:,2); % oxygen volume percent rear
CO2_rear_fds_200=RSE_200data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_200=RSE_200data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_200=RSE_200data(:,5); % fuel volume percent rear
H2O_rear_fds_200=RSE_200data(:,6); % water vapor volume percent rear

O2_front_fds_200=RSE_200data(:,7); % oxygen volume percent front
CO2_front_fds_200=RSE_200data(:,8); % carbon dioxide volume percent front
CO_front_fds_200=RSE_200data(:,9); % carbon monoxide volume percent front
Methane_front_fds_200=RSE_200data(:,10); % fuel volume percent front
H2O_front_fds_200=RSE_200data(:,11); % water vapor volume percent front

temp_rearA_fds_200=RSE_200data(:,12); % temperature rear sample A
temp_rearB_fds_200=RSE_200data(:,13); % temperature rear sample B

temp_frontA_fds_200=RSE_200data(:,14); % temperature front sample A
temp_frontB_fds_200=RSE_200data(:,15); % temperature front sample B
% -------------
% 300 kW
% -------------
RSE_300_struct=importdata('NIST_RSE_1994_300_devc.csv'); %FDS data file
RSE_300_header=RSE_300_struct.textdata; % headers
RSE_300data=RSE_300_struct.data; % data

time300=RSE_300data(:,1); % time

O2_rear_fds_300=RSE_300data(:,2); % oxygen volume percent rear
CO2_rear_fds_300=RSE_300data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_300=RSE_300data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_300=RSE_300data(:,5); % fuel volume percent rear
H2O_rear_fds_300=RSE_300data(:,6); % water vapor volume percent rear

O2_front_fds_300=RSE_300data(:,7); % oxygen volume percent front
CO2_front_fds_300=RSE_300data(:,8); % carbon dioxide volume percent front
CO_front_fds_300=RSE_300data(:,9); % carbon monoxide volume percent front
Methane_front_fds_300=RSE_300data(:,10); % fuel volume percent front
H2O_front_fds_300=RSE_300data(:,11); % water vapor volume percent front

temp_rearA_fds_300=RSE_300data(:,12); % temperature rear sample A
temp_rearB_fds_300=RSE_300data(:,13); % temperature rear sample B

temp_frontA_fds_300=RSE_300data(:,14); % temperature front sample A
temp_frontB_fds_300=RSE_300data(:,15); % temperature front sample B
% -------------
% 400 kW
% -------------
RSE_400_struct=importdata('NIST_RSE_1994_400_devc.csv'); %FDS data file
RSE_400_header=RSE_400_struct.textdata; % headers
RSE_400data=RSE_400_struct.data; % data

time400=RSE_400data(:,1); % time

O2_rear_fds_400=RSE_400data(:,2); % oxygen volume percent rear
CO2_rear_fds_400=RSE_400data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_400=RSE_400data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_400=RSE_400data(:,5); % fuel volume percent rear
H2O_rear_fds_400=RSE_400data(:,6); % water vapor volume percent rear

O2_front_fds_400=RSE_400data(:,7); % oxygen volume percent front
CO2_front_fds_400=RSE_400data(:,8); % carbon dioxide volume percent front
CO_front_fds_400=RSE_400data(:,9); % carbon monoxide volume percent front
Methane_front_fds_400=RSE_400data(:,10); % fuel volume percent front
H2O_front_fds_400=RSE_400data(:,11); % water vapor volume percent front

temp_rearA_fds_400=RSE_400data(:,12); % temperature rear sample A
temp_rearB_fds_400=RSE_400data(:,13); % temperature rear sample B

temp_frontA_fds_400=RSE_400data(:,14); % temperature front sample A
temp_frontB_fds_400=RSE_400data(:,15); % temperature front sample B
% -------------
% 500 kW
% -------------
RSE_500_struct=importdata('NIST_RSE_1994_500_devc.csv'); %FDS data file
RSE_500_header=RSE_500_struct.textdata; % headers
RSE_500data=RSE_500_struct.data; % data

time500=RSE_500data(:,1); % time

O2_rear_fds_500=RSE_500data(:,2); % oxygen volume percent rear
CO2_rear_fds_500=RSE_500data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_500=RSE_500data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_500=RSE_500data(:,5); % fuel volume percent rear
H2O_rear_fds_500=RSE_500data(:,6); % water vapor volume percent rear

O2_front_fds_500=RSE_500data(:,7); % oxygen volume percent front
CO2_front_fds_500=RSE_500data(:,8); % carbon dioxide volume percent front
CO_front_fds_500=RSE_500data(:,9); % carbon monoxide volume percent front
Methane_front_fds_500=RSE_500data(:,10); % fuel volume percent front
H2O_front_fds_500=RSE_500data(:,11); % water vapor volume percent front

temp_rearA_fds_500=RSE_500data(:,12); % temperature rear sample A
temp_rearB_fds_500=RSE_500data(:,13); % temperature rear sample B

temp_frontA_fds_500=RSE_500data(:,14); % temperature front sample A
temp_frontB_fds_500=RSE_500data(:,15); % temperature front sample B
% -------------
% 600 kW
% -------------
RSE_600_struct=importdata('NIST_RSE_1994_600_devc.csv'); %FDS data file
RSE_600_header=RSE_600_struct.textdata; % headers
RSE_600data=RSE_600_struct.data; % data

time600=RSE_600data(:,1); % time

O2_rear_fds_600=RSE_600data(:,2); % oxygen volume percent rear
CO2_rear_fds_600=RSE_600data(:,3); % carbon dioxide volume percent rear
CO_rear_fds_600=RSE_600data(:,4); % carbon monoxide volume percent rear
Methane_rear_fds_600=RSE_600data(:,5); % fuel volume percent rear
H2O_rear_fds_600=RSE_600data(:,6); % water vapor volume percent rear

O2_front_fds_600=RSE_600data(:,7); % oxygen volume percent front
CO2_front_fds_600=RSE_600data(:,8); % carbon dioxide volume percent front
CO_front_fds_600=RSE_600data(:,9); % carbon monoxide volume percent front
Methane_front_fds_600=RSE_600data(:,10); % fuel volume percent front
H2O_front_fds_600=RSE_600data(:,11); % water vapor volume percent front

temp_rearA_fds_600=RSE_600data(:,12); % temperature rear sample A
temp_rearB_fds_600=RSE_600data(:,13); % temperature rear sample B

temp_frontA_fds_600=RSE_600data(:,14); % temperature front sample A
temp_frontB_fds_600=RSE_600data(:,15); % temperature front sample B

%------------------------------------------------
% Building FDS data file with average species/temp
% Function of HRR
%------------------------------------------------
start_row=12; % row to begin averaging

HRR= [50 75 100 150 200 300 400 500 600];

%-------------
% Mean O2 rear
%-------------
meanO2_rear(1,1)=mean(O2_rear_fds_50(start_row:end));
meanO2_rear(2,1)=mean(O2_rear_fds_75(start_row:end));
meanO2_rear(3,1)=mean(O2_rear_fds_100(start_row:end));
meanO2_rear(4,1)=mean(O2_rear_fds_150(start_row:end));
meanO2_rear(5,1)=mean(O2_rear_fds_200(start_row:end));
meanO2_rear(6,1)=mean(O2_rear_fds_300(start_row:end));
meanO2_rear(7,1)=mean(O2_rear_fds_400(start_row:end));
meanO2_rear(8,1)=mean(O2_rear_fds_500(start_row:end));
meanO2_rear(9,1)=mean(O2_rear_fds_600(start_row:end));
%-------------
% Mean O2 front
%-------------
meanO2_front(1,1)=mean(O2_front_fds_50(start_row:end));
meanO2_front(2,1)=mean(O2_front_fds_75(start_row:end));
meanO2_front(3,1)=mean(O2_front_fds_100(start_row:end));
meanO2_front(4,1)=mean(O2_front_fds_150(start_row:end));
meanO2_front(5,1)=mean(O2_front_fds_200(start_row:end));
meanO2_front(6,1)=mean(O2_front_fds_300(start_row:end));
meanO2_front(7,1)=mean(O2_front_fds_400(start_row:end));
meanO2_front(8,1)=mean(O2_front_fds_500(start_row:end));
meanO2_front(9,1)=mean(O2_front_fds_600(start_row:end));
%-------------
% Mean CO2 rear
%-------------
meanCO2_rear(1,1)=mean(CO2_rear_fds_50(start_row:end));
meanCO2_rear(2,1)=mean(CO2_rear_fds_75(start_row:end));
meanCO2_rear(3,1)=mean(CO2_rear_fds_100(start_row:end));
meanCO2_rear(4,1)=mean(CO2_rear_fds_150(start_row:end));
meanCO2_rear(5,1)=mean(CO2_rear_fds_200(start_row:end));
meanCO2_rear(6,1)=mean(CO2_rear_fds_300(start_row:end));
meanCO2_rear(7,1)=mean(CO2_rear_fds_400(start_row:end));
meanCO2_rear(8,1)=mean(CO2_rear_fds_500(start_row:end));
meanCO2_rear(9,1)=mean(CO2_rear_fds_600(start_row:end));
%-------------
% Mean CO2 front
%-------------
meanCO2_front(1,1)=mean(CO2_front_fds_50(start_row:end));
meanCO2_front(2,1)=mean(CO2_front_fds_75(start_row:end));
meanCO2_front(3,1)=mean(CO2_front_fds_100(start_row:end));
meanCO2_front(4,1)=mean(CO2_front_fds_150(start_row:end));
meanCO2_front(5,1)=mean(CO2_front_fds_200(start_row:end));
meanCO2_front(6,1)=mean(CO2_front_fds_300(start_row:end));
meanCO2_front(7,1)=mean(CO2_front_fds_400(start_row:end));
meanCO2_front(8,1)=mean(CO2_front_fds_500(start_row:end));
meanCO2_front(9,1)=mean(CO2_front_fds_600(start_row:end));
%-------------
% Mean CO rear
%-------------
meanCO_rear(1,1)=mean(CO_rear_fds_50(start_row:end));
meanCO_rear(2,1)=mean(CO_rear_fds_75(start_row:end));
meanCO_rear(3,1)=mean(CO_rear_fds_100(start_row:end));
meanCO_rear(4,1)=mean(CO_rear_fds_150(start_row:end));
meanCO_rear(5,1)=mean(CO_rear_fds_200(start_row:end));
meanCO_rear(6,1)=mean(CO_rear_fds_300(start_row:end));
meanCO_rear(7,1)=mean(CO_rear_fds_400(start_row:end));
meanCO_rear(8,1)=mean(CO_rear_fds_500(start_row:end));
meanCO_rear(9,1)=mean(CO_rear_fds_600(start_row:end));
%-------------
% Mean CO front
%-------------
meanCO_front(1,1)=mean(CO_front_fds_50(start_row:end));
meanCO_front(2,1)=mean(CO_front_fds_75(start_row:end));
meanCO_front(3,1)=mean(CO_front_fds_100(start_row:end));
meanCO_front(4,1)=mean(CO_front_fds_150(start_row:end));
meanCO_front(5,1)=mean(CO_front_fds_200(start_row:end));
meanCO_front(6,1)=mean(CO_front_fds_300(start_row:end));
meanCO_front(7,1)=mean(CO_front_fds_400(start_row:end));
meanCO_front(8,1)=mean(CO_front_fds_500(start_row:end));
meanCO_front(9,1)=mean(CO_front_fds_600(start_row:end));
%-------------
% Mean Methane rear
%-------------
meanMethane_rear(1,1)=mean(Methane_rear_fds_50(start_row:end));
meanMethane_rear(2,1)=mean(Methane_rear_fds_75(start_row:end));
meanMethane_rear(3,1)=mean(Methane_rear_fds_100(start_row:end));
meanMethane_rear(4,1)=mean(Methane_rear_fds_150(start_row:end));
meanMethane_rear(5,1)=mean(Methane_rear_fds_200(start_row:end));
meanMethane_rear(6,1)=mean(Methane_rear_fds_300(start_row:end));
meanMethane_rear(7,1)=mean(Methane_rear_fds_400(start_row:end));
meanMethane_rear(8,1)=mean(Methane_rear_fds_500(start_row:end));
meanMethane_rear(9,1)=mean(Methane_rear_fds_600(start_row:end));
%-------------
% Mean Methane front
%-------------
meanMethane_front(1,1)=mean(Methane_front_fds_50(start_row:end));
meanMethane_front(2,1)=mean(Methane_front_fds_75(start_row:end));
meanMethane_front(3,1)=mean(Methane_front_fds_100(start_row:end));
meanMethane_front(4,1)=mean(Methane_front_fds_150(start_row:end));
meanMethane_front(5,1)=mean(Methane_front_fds_200(start_row:end));
meanMethane_front(6,1)=mean(Methane_front_fds_300(start_row:end));
meanMethane_front(7,1)=mean(Methane_front_fds_400(start_row:end));
meanMethane_front(8,1)=mean(Methane_front_fds_500(start_row:end));
meanMethane_front(9,1)=mean(Methane_front_fds_600(start_row:end));
%-------------
% Mean TempA rear
%-------------
meanTRA(1,1)=mean(temp_rearA_fds_50(start_row:end));
meanTRA(2,1)=mean(temp_rearA_fds_75(start_row:end));
meanTRA(3,1)=mean(temp_rearA_fds_100(start_row:end));
meanTRA(4,1)=mean(temp_rearA_fds_150(start_row:end));
meanTRA(5,1)=mean(temp_rearA_fds_200(start_row:end));
meanTRA(6,1)=mean(temp_rearA_fds_300(start_row:end));
meanTRA(7,1)=mean(temp_rearA_fds_400(start_row:end));
meanTRA(8,1)=mean(temp_rearA_fds_500(start_row:end));
meanTRA(9,1)=mean(temp_rearA_fds_600(start_row:end));
%-------------
% Mean TempA front
%-------------
meanTFA(1,1)=mean(temp_frontA_fds_50(start_row:end));
meanTFA(2,1)=mean(temp_frontA_fds_75(start_row:end));
meanTFA(3,1)=mean(temp_frontA_fds_100(start_row:end));
meanTFA(4,1)=mean(temp_frontA_fds_150(start_row:end));
meanTFA(5,1)=mean(temp_frontA_fds_200(start_row:end));
meanTFA(6,1)=mean(temp_frontA_fds_300(start_row:end));
meanTFA(7,1)=mean(temp_frontA_fds_400(start_row:end));
meanTFA(8,1)=mean(temp_frontA_fds_500(start_row:end));
meanTFA(9,1)=mean(temp_frontA_fds_600(start_row:end));
%-------------
% Mean TempB rear
%-------------
meanTRBB(1,1)=mean(temp_rearB_fds_50(start_row:end));
meanTRBB(2,1)=mean(temp_rearB_fds_75(start_row:end));
meanTRBB(3,1)=mean(temp_rearB_fds_100(start_row:end));
meanTRBB(4,1)=mean(temp_rearB_fds_150(start_row:end));
meanTRBB(5,1)=mean(temp_rearB_fds_200(start_row:end));
meanTRBB(6,1)=mean(temp_rearB_fds_300(start_row:end));
meanTRBB(7,1)=mean(temp_rearB_fds_400(start_row:end));
meanTRBB(8,1)=mean(temp_rearB_fds_500(start_row:end));
meanTRBB(9,1)=mean(temp_rearB_fds_600(start_row:end));
%-------------
% Mean TempB front
%-------------
meanTFBB(1,1)=mean(temp_frontB_fds_50(start_row:end));
meanTFBB(2,1)=mean(temp_frontB_fds_75(start_row:end));
meanTFBB(3,1)=mean(temp_frontB_fds_100(start_row:end));
meanTFBB(4,1)=mean(temp_frontB_fds_150(start_row:end));
meanTFBB(5,1)=mean(temp_frontB_fds_200(start_row:end));
meanTFBB(6,1)=mean(temp_frontB_fds_300(start_row:end));
meanTFBB(7,1)=mean(temp_frontB_fds_400(start_row:end));
meanTFBB(8,1)=mean(temp_frontB_fds_500(start_row:end));
meanTFBB(9,1)=mean(temp_frontB_fds_600(start_row:end));

%------------------------------------------------
% FDS Data CSV File
%------------------------------------------------

RSE_Results(:,1)=HRR;
RSE_Results(:,2)=meanO2_rear;
RSE_Results(:,3)=meanCO2_rear;
RSE_Results(:,4)=meanCO_rear;
RSE_Results(:,5)=meanMethane_rear;
RSE_Results(:,6)=meanO2_front;
RSE_Results(:,7)=meanCO2_front;
RSE_Results(:,8)=meanCO_front;
RSE_Results(:,9)=meanMethane_front;
RSE_Results(:,10)=meanTRA;
RSE_Results(:,11)=meanTRBB;
RSE_Results(:,12)=meanTFA;
RSE_Results(:,13)=meanTFBB;

header1 = {'HRR','O2Rear_FDS','CO2Rear_FDS','CORear_FDS','UHRear_FDS','O2Front_FDS','CO2Front_FDS',...
    'COFront_FDS','UHFront_FDS','TRSampA_FDS','TRSampBB_FDS','TFSampA_FDS','TFSampBB_FDS'};
filename1 = '../../Validation/NIST_RSE_1994//FDS_Output_Files/NIST_RSE_1994_FDS.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',header1{:});
for i=1:9
    fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',RSE_Results(i,:));
end
fclose(fid);


