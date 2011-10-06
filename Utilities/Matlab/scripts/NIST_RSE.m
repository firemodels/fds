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

addpath('../../Validation/NIST_RSE_1994/FDS_Output_Files/')

%------------------------------------------------
% Read in FDS output files
% 9 FDS HRRs
%------------------------------------------------

%-------------
% Setup FDS Device Headers
%-------------
Header = {'Time','"O2Rear_FDS"','"CO2Rear_FDS"','"CORear_FDS"','"UHRear_FDS"',...
    '"H2ORear_FDS"','"O2Front_FDS"','"CO2Front_FDS"','"COFront_FDS"','"UHFront_FDS"',...
    '"H2OFront_FDS"','"TRSampA_FDS"','"TRSampBB_FDS"','"TFSampA_FDS"','"TFSampBB_FDS"'};

%-------------
% 50 kW
%-------------
RSE_50_struct=importdata('NIST_RSE_1994_50_devc.csv'); %FDS data file
RSE_50_header=RSE_50_struct.textdata; % headers
RSE_50data=RSE_50_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_50_header(2,i)) 
        time50=RSE_50data(:,i); % time
    elseif strcmp(Header(2),RSE_50_header(2,i))
        O2_rear_fds_50=RSE_50data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_50_header(2,i))
        CO2_rear_fds_50=RSE_50data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_50_header(2,i))
        CO_rear_fds_50=RSE_50data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_50_header(2,i))
        Methane_rear_fds_50=RSE_50data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_50_header(2,i))
        H2O_rear_fds_50=RSE_50data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_50_header(2,i))
        O2_front_fds_50=RSE_50data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_50_header(2,i))
        CO2_front_fds_50=RSE_50data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_50_header(2,i))
        CO_front_fds_50=RSE_50data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_50_header(2,i))
        Methane_front_fds_50=RSE_50data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_50_header(2,i))
        H2O_front_fds_50=RSE_50data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_50_header(2,i))
        temp_rearA_fds_50=RSE_50data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_50_header(2,i))
        temp_rearB_fds_50=RSE_50data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_50_header(2,i))
        temp_frontA_fds_50=RSE_50data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_50_header(2,i))
        temp_frontB_fds_50=RSE_50data(:,i); % temperature front sample
    end
end
% -------------
% 75 kW
% -------------
RSE_75_struct=importdata('NIST_RSE_1994_75_devc.csv'); %FDS data file
RSE_75_header=RSE_75_struct.textdata; % headers
RSE_75data=RSE_75_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_75_header(2,i)) 
        time75=RSE_75data(:,i); % time
    elseif strcmp(Header(2),RSE_75_header(2,i))
        O2_rear_fds_75=RSE_75data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_75_header(2,i))
        CO2_rear_fds_75=RSE_75data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_75_header(2,i))
        CO_rear_fds_75=RSE_75data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_75_header(2,i))
        Methane_rear_fds_75=RSE_75data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_75_header(2,i))
        H2O_rear_fds_75=RSE_75data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_75_header(2,i))
        O2_front_fds_75=RSE_75data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_75_header(2,i))
        CO2_front_fds_75=RSE_75data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_75_header(2,i))
        CO_front_fds_75=RSE_75data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_75_header(2,i))
        Methane_front_fds_75=RSE_75data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_75_header(2,i))
        H2O_front_fds_75=RSE_75data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_75_header(2,i))
        temp_rearA_fds_75=RSE_75data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_75_header(2,i))
        temp_rearB_fds_75=RSE_75data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_75_header(2,i))
        temp_frontA_fds_75=RSE_75data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_75_header(2,i))
        temp_frontB_fds_75=RSE_75data(:,i); % temperature front sample
    end
end
% -------------
% 100 kW
% -------------
RSE_100_struct=importdata('NIST_RSE_1994_100_devc.csv'); %FDS data file
RSE_100_header=RSE_100_struct.textdata; % headers
RSE_100data=RSE_100_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_100_header(2,i)) 
        time100=RSE_100data(:,i); % time
    elseif strcmp(Header(2),RSE_100_header(2,i))
        O2_rear_fds_100=RSE_100data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_100_header(2,i))
        CO2_rear_fds_100=RSE_100data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_100_header(2,i))
        CO_rear_fds_100=RSE_100data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_100_header(2,i))
        Methane_rear_fds_100=RSE_100data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_100_header(2,i))
        H2O_rear_fds_100=RSE_100data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_100_header(2,i))
        O2_front_fds_100=RSE_100data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_100_header(2,i))
        CO2_front_fds_100=RSE_100data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_100_header(2,i))
        CO_front_fds_100=RSE_100data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_100_header(2,i))
        Methane_front_fds_100=RSE_100data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_100_header(2,i))
        H2O_front_fds_100=RSE_100data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_100_header(2,i))
        temp_rearA_fds_100=RSE_100data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_100_header(2,i))
        temp_rearB_fds_100=RSE_100data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_100_header(2,i))
        temp_frontA_fds_100=RSE_100data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_100_header(2,i))
        temp_frontB_fds_100=RSE_100data(:,i); % temperature front sample
    end
end
% -------------
% 150 kW
% -------------
RSE_150_struct=importdata('NIST_RSE_1994_150_devc.csv'); %FDS data file
RSE_150_header=RSE_150_struct.textdata; % headers
RSE_150data=RSE_150_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_150_header(2,i)) 
        time150=RSE_150data(:,i); % time
    elseif strcmp(Header(2),RSE_150_header(2,i))
        O2_rear_fds_150=RSE_150data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_150_header(2,i))
        CO2_rear_fds_150=RSE_150data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_150_header(2,i))
        CO_rear_fds_150=RSE_150data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_150_header(2,i))
        Methane_rear_fds_150=RSE_150data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_150_header(2,i))
        H2O_rear_fds_150=RSE_150data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_150_header(2,i))
        O2_front_fds_150=RSE_150data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_150_header(2,i))
        CO2_front_fds_150=RSE_150data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_150_header(2,i))
        CO_front_fds_150=RSE_150data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_150_header(2,i))
        Methane_front_fds_150=RSE_150data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_150_header(2,i))
        H2O_front_fds_150=RSE_150data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_150_header(2,i))
        temp_rearA_fds_150=RSE_150data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_150_header(2,i))
        temp_rearB_fds_150=RSE_150data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_150_header(2,i))
        temp_frontA_fds_150=RSE_150data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_150_header(2,i))
        temp_frontB_fds_150=RSE_150data(:,i); % temperature front sample
    end
end
% -------------
% 200 kW
% -------------
RSE_200_struct=importdata('NIST_RSE_1994_200_devc.csv'); %FDS data file
RSE_200_header=RSE_200_struct.textdata; % headers
RSE_200data=RSE_200_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_200_header(2,i)) 
        time200=RSE_200data(:,i); % time
    elseif strcmp(Header(2),RSE_200_header(2,i))
        O2_rear_fds_200=RSE_200data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_200_header(2,i))
        CO2_rear_fds_200=RSE_200data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_200_header(2,i))
        CO_rear_fds_200=RSE_200data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_200_header(2,i))
        Methane_rear_fds_200=RSE_200data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_200_header(2,i))
        H2O_rear_fds_200=RSE_200data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_200_header(2,i))
        O2_front_fds_200=RSE_200data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_200_header(2,i))
        CO2_front_fds_200=RSE_200data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_200_header(2,i))
        CO_front_fds_200=RSE_200data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_200_header(2,i))
        Methane_front_fds_200=RSE_200data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_200_header(2,i))
        H2O_front_fds_200=RSE_200data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_200_header(2,i))
        temp_rearA_fds_200=RSE_200data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_200_header(2,i))
        temp_rearB_fds_200=RSE_200data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_200_header(2,i))
        temp_frontA_fds_200=RSE_200data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_200_header(2,i))
        temp_frontB_fds_200=RSE_200data(:,i); % temperature front sample
    end
end
% -------------
% 300 kW
% -------------
RSE_300_struct=importdata('NIST_RSE_1994_300_devc.csv'); %FDS data file
RSE_300_header=RSE_300_struct.textdata; % headers
RSE_300data=RSE_300_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_300_header(2,i)) 
        time300=RSE_300data(:,i); % time
    elseif strcmp(Header(2),RSE_300_header(2,i))
        O2_rear_fds_300=RSE_300data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_300_header(2,i))
        CO2_rear_fds_300=RSE_300data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_300_header(2,i))
        CO_rear_fds_300=RSE_300data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_300_header(2,i))
        Methane_rear_fds_300=RSE_300data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_300_header(2,i))
        H2O_rear_fds_300=RSE_300data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_300_header(2,i))
        O2_front_fds_300=RSE_300data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_300_header(2,i))
        CO2_front_fds_300=RSE_300data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_300_header(2,i))
        CO_front_fds_300=RSE_300data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_300_header(2,i))
        Methane_front_fds_300=RSE_300data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_300_header(2,i))
        H2O_front_fds_300=RSE_300data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_300_header(2,i))
        temp_rearA_fds_300=RSE_300data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_300_header(2,i))
        temp_rearB_fds_300=RSE_300data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_300_header(2,i))
        temp_frontA_fds_300=RSE_300data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_300_header(2,i))
        temp_frontB_fds_300=RSE_300data(:,i); % temperature front sample
    end
end
% -------------
% 400 kW
% -------------
RSE_400_struct=importdata('NIST_RSE_1994_400_devc.csv'); %FDS data file
RSE_400_header=RSE_400_struct.textdata; % headers
RSE_400data=RSE_400_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_400_header(2,i)) 
        time400=RSE_400data(:,i); % time
    elseif strcmp(Header(2),RSE_400_header(2,i))
        O2_rear_fds_400=RSE_400data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_400_header(2,i))
        CO2_rear_fds_400=RSE_400data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_400_header(2,i))
        CO_rear_fds_400=RSE_400data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_400_header(2,i))
        Methane_rear_fds_400=RSE_400data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_400_header(2,i))
        H2O_rear_fds_400=RSE_400data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_400_header(2,i))
        O2_front_fds_400=RSE_400data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_400_header(2,i))
        CO2_front_fds_400=RSE_400data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_400_header(2,i))
        CO_front_fds_400=RSE_400data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_400_header(2,i))
        Methane_front_fds_400=RSE_400data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_400_header(2,i))
        H2O_front_fds_400=RSE_400data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_400_header(2,i))
        temp_rearA_fds_400=RSE_400data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_400_header(2,i))
        temp_rearB_fds_400=RSE_400data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_400_header(2,i))
        temp_frontA_fds_400=RSE_400data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_400_header(2,i))
        temp_frontB_fds_400=RSE_400data(:,i); % temperature front sample
    end
end
% -------------
% 500 kW
% -------------
RSE_500_struct=importdata('NIST_RSE_1994_500_devc.csv'); %FDS data file
RSE_500_header=RSE_500_struct.textdata; % headers
RSE_500data=RSE_500_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_500_header(2,i)) 
        time500=RSE_500data(:,i); % time
    elseif strcmp(Header(2),RSE_500_header(2,i))
        O2_rear_fds_500=RSE_500data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_500_header(2,i))
        CO2_rear_fds_500=RSE_500data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_500_header(2,i))
        CO_rear_fds_500=RSE_500data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_500_header(2,i))
        Methane_rear_fds_500=RSE_500data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_500_header(2,i))
        H2O_rear_fds_500=RSE_500data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_500_header(2,i))
        O2_front_fds_500=RSE_500data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_500_header(2,i))
        CO2_front_fds_500=RSE_500data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_500_header(2,i))
        CO_front_fds_500=RSE_500data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_500_header(2,i))
        Methane_front_fds_500=RSE_500data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_500_header(2,i))
        H2O_front_fds_500=RSE_500data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_500_header(2,i))
        temp_rearA_fds_500=RSE_500data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_500_header(2,i))
        temp_rearB_fds_500=RSE_500data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_500_header(2,i))
        temp_frontA_fds_500=RSE_500data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_500_header(2,i))
        temp_frontB_fds_500=RSE_500data(:,i); % temperature front sample
    end
end
% -------------
% 600 kW
% -------------
RSE_600_struct=importdata('NIST_RSE_1994_600_devc.csv'); %FDS data file
RSE_600_header=RSE_600_struct.textdata; % headers
RSE_600data=RSE_600_struct.data; % data

for i=1:length(Header)
    if strcmp(Header(1),RSE_600_header(2,i)) 
        time600=RSE_600data(:,i); % time
    elseif strcmp(Header(2),RSE_600_header(2,i))
        O2_rear_fds_600=RSE_600data(:,i); % oxygen volume percent rear
    elseif strcmp(Header(3),RSE_600_header(2,i))
        CO2_rear_fds_600=RSE_600data(:,i); % carbon dioxide volume percent rear
    elseif strcmp(Header(4),RSE_600_header(2,i))
        CO_rear_fds_600=RSE_600data(:,i); % carbon monoxide volume percent rear
    elseif strcmp(Header(5),RSE_600_header(2,i))
        Methane_rear_fds_600=RSE_600data(:,i); % fuel volume percent rear
    elseif strcmp(Header(6),RSE_600_header(2,i))
        H2O_rear_fds_600=RSE_600data(:,i); % water vapor volume percent rear
    elseif strcmp(Header(7),RSE_600_header(2,i))
        O2_front_fds_600=RSE_600data(:,i); % oxygen volume percent front
    elseif strcmp(Header(8),RSE_600_header(2,i))
        CO2_front_fds_600=RSE_600data(:,i); % carbon dioxide volume percent front
    elseif strcmp(Header(9),RSE_600_header(2,i))
        CO_front_fds_600=RSE_600data(:,i); % carbon monoxide volume percent front
    elseif strcmp(Header(10),RSE_600_header(2,i))
        Methane_front_fds_600=RSE_600data(:,i); % fuel volume percent front
    elseif strcmp(Header(11),RSE_600_header(2,i))
        H2O_front_fds_600=RSE_600data(:,i); % water vapor volume percent front
    elseif strcmp(Header(12),RSE_600_header(2,i))
        temp_rearA_fds_600=RSE_600data(:,i); % temperature rear sample A
    elseif strcmp(Header(13),RSE_600_header(2,i))
        temp_rearB_fds_600=RSE_600data(:,i); % temperature rear sample B
    elseif strcmp(Header(14),RSE_600_header(2,i))
        temp_frontA_fds_600=RSE_600data(:,i); % temperature front sample A
    elseif strcmp(Header(15),RSE_600_header(2,i))
        temp_frontB_fds_600=RSE_600data(:,i); % temperature front sample
    end
end

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
filename1 = '../../Validation/NIST_RSE_1994/FDS_Output_Files/NIST_RSE_1994_FDS.csv';
fid = fopen(filename1,'wt');
fprintf(fid,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',header1{:});
for i=1:9
    fprintf(fid,'%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',RSE_Results(i,:));
end
fclose(fid);
