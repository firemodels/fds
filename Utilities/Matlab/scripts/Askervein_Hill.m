% N. Crump  6/11/2021 - 8/5/2021
%Matlabs script to transpose the Askervein hill case files into a form
%usable by the FDS validation script.

% Script Process:
%-Read in the decv.csv file from the out directory
%-Transpose the matrix such that each device has it's own row, not column
%-Identify the ASW and ANE (Line A SouthWest -> NorthEast)Devices
%-Insert the Distance from Hiltop column appropriately
%-Sort Rows in Acending order (-100->40)
%-output new devc.csv file

% Navigation
% to jump to each file export in order, use crtl-f and search for 'NewName ='
%-10m Post lines A, AA, B
     %ASW-ANE:   SPEED DIRECTION UPWASH SIGM_SPEED SIGMw
     %BNW-BSE:   SPEED
     %AASW-AANE: SPEED
%-Mean flow runs MF and TU
     %TU:   SPEED DIRECTION UPWASH SIGM_SPEED SIGMw
     %MF:   SPEED
%-50m Profile at HT 
   %     Data:   SPEED, SIGM_SPEED
%-FRG 17m Profile at CP' 
     %TU:   SPEED DIRECTION UPWASH SIGM_SPEED
     %MF:   SPEED
%-UK 30m Profile at ASW60 
   %     Data:   SPEED DIRECTION UPWASH SIGM_SPEED

%Questions Checklist:
% -Why did blaze add "..." to the device titles? Is it my fault? it makes
% reading in the column headers very strange -> '"..."' but the code will
% no longer work without them. To check this, >> A.colheaders
 
%Debugging notes:
% -The most common error will be dimensions of arrays cannot be concacted
% this is most often caused by a failure to read in the the script. 
%   Check if Time_Stop is set to the correct value
%   Check if the file has been saved in column veiwing format
%   Check if the colheaders have '".."' or '..' with >> A.colheaders
%   Check values with >> whos DIST DIRECTION UPWASH AVG_SPEED SIGM_SPEED SIGMw


% outdir = '../../../../out/Askervein_Hill/';
outdir = '../../../out/Askervein_Hill/';

%import data from file
filename = 'Askervein_TU03A_16m_devc.csv';
A = importdata([outdir,filename], ",", 2);

% isolates unique part of device file name for easeir resolution comparisons
 ID = strrep(strrep(filename,'Askervein_TU03A_',''),'_devc.csv','');

%Time_Stop allows the user to change what time they want from the table
%It's currently set up to automatically take the latest avarege
% TimeStop = 3600
Times = A.data(:,(strcmp(A.colheaders,'Time')));
TimeStop = max(Times);
TimeStop1 = TimeStop -1;
TimeStop2 = TimeStop +1;
slot = find( A.data(:,(strcmp(A.colheaders,'Time'))) > TimeStop1 & ...
             A.data(:,(strcmp(A.colheaders,'Time'))) < TimeStop2 );

% A Line 10m height
% ----------------------------------------------------------------
NewName =strcat('ASW-ANE_TU03A_',ID,'.csv');
% Create Distance from hiltop column. 
DIST = ["(m)", "DIST", -85, -60, -50, -35, -20, -10, 0, 10, 20, 40]';
         
%data at TIME_STOP (slot variable), this array holds the avarage SPEED data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT 10 m t"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW UK 30 m tower 10m"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ANE40"')));
SPEED = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

%data at TIME_STOP (slot variable), this array holds the SIGM SPEED data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT 10 m t"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW UK 30 m tower 10m"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ANE40"')));
% SIGM_SPEED is not used to calculate other values, so header is added here
SIGM_SPEED = ["(m/s)","SIGM_SPEED",ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,]';

%data at TIME_STOP (slot variable), this array holds the SIGMw data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"SIGMw HT 10 m t"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW UK 30 m tower 10m"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ANE40"')));
SIGMw = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];
INTw = SIGMw./SPEED;
SIGMw = ["(m/s)";"SIGMw";SIGMw'];

%data at TIME_STOP (slot variable), this array holds the AVGu data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"AVGu HT 10 m t"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW UK 30 m tower 10m"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ANE40"')));
AVGu = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

%data at TIME_STOP (slot variable), this array holds the AVGv data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"AVGv HT 10 m t"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW UK 30 m tower 10m"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ANE40"')));
AVGv = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

%data at TIME_STOP (slot variable), this array holds the AVGw data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"AVGw HT 10 m t"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW UK 30 m tower 10m"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ANE40"')));
AVGw = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

% calculate secondary variables and add headers to primary variable columns
HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
AVG_SPEED = ["(m/s)";"SPEED";SPEED'];

% Create final table
FINAL = [DIST DIRECTION UPWASH AVG_SPEED SIGM_SPEED SIGMw];
% Write final table to CSV file in OUT directory
writematrix(FINAL,[outdir, NewName])

% AA Line
% ----------------------------------------------------------------
%SPEED
NewName =strcat('AASW-AANE_TU03A_',ID,'.csv');
%NewName ='AASW-AANE_TU03A_10m.csv';
DIST = ["(m)","DIST",-90,-80,-70,-60,-50,-40,-30,-20,-10,0,0,10,20,30,40,50,60]';
CP_BSE40 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP BSE 40"')));
AACP0 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP UK mf"')));
% AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW10 t"')));
% AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW30 t"')));
% AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW50 t"')));
AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW10 mf"')));
AASW20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW20"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW30 mf"')));
AASW40 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW40"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW50 mf"')));
AASW60 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW60"')));
AASW70 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW70"')));
AASW80 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW80"')));
AASW90 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW90"')));
SPEED = [AASW90,AASW80,AASW70,AASW60,AASW50,AASW40,AASW30,AASW20,AASW10];
AANE10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AANE10"')));
AANE20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AANE20"')));
AANE30 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AANE30"')));
AANE40 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AANE40"')));
AANE50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AANE50"')));
AANE60 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AANE60"')));
SPEED = [SPEED,CP_BSE40,AACP0,AANE10,AANE20,AANE30,AANE40,AANE50,AANE60,];
AVG_SPEED = ["(m/s)";"SPEED";SPEED'];

FINAL = [DIST AVG_SPEED];
writematrix(FINAL,[outdir, NewName])

% B Line
% ----------------------------------------------------------------
%SPEED
NewName =strcat('BNW-BSE_TU03A_',ID,'.csv');
% NewName ='BNW-BSE_TU03A_10m.csv';
DIST = ["(m)","DIST",-20,-10,0,10,20,30,40,40,50,60,70,80,90,100,110,130]';
HT_0 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT 10 m mf"')));
BNW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BNW10"')));
BNW20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BNW20"')));
BSE10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE10"')));
BSE20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE20"')));
BSE30 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE30"')));
BSE40 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE40"')));
BSE50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE50"')));
BSE60 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE60"')));
BSE70 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE70"')));
BSE80 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE80"')));
BSE90 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE90"')));
BSE100 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE100"')));
BSE110 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE110"')));
BSE130 = A.data(slot,find(strcmp(A.colheaders,'"SPEED BSE130"')));
AACP0 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP UK mf"')));

SPEED = [BNW20,BNW10,HT_0,BSE10,BSE20,BSE30,BSE40,AACP0,BSE50,BSE60,BSE70,BSE80,BSE90,BSE100,BSE110,BSE130];
AVG_SPEED = ["(m/s)";"SPEED";SPEED'];
FINAL = [DIST AVG_SPEED];
writematrix(FINAL,[outdir, NewName])

% AA Short Line mf and tu Mean Flow
% ----------------------------------------------------------------
%AA MF
NewName =strcat('AASW-AANE_mf_TU03A_',ID,'.csv');
% NewName ='AASW-AANE_mf_TU03A.csv';
DIST = ["(m)", "DIST", -50, -30, -10, 0]';
AASW0 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG mf"')));
AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW10 mf"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW30 mf"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW50 mf"')));
SPEED = [AASW0,AASW10,AASW30,AASW50];

AASW0 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG mf"')));
AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED AASW10 mf"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED AASW30 mf"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED AASW50 mf"')));
SIGM_SPEED = ["(m/s)","SIGM_SPEED",AASW0,AASW10,AASW30,AASW50]';

AASW0 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG mf"')));
AASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu AASW10 mf"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"AVGu AASW30 mf"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGu AASW50 mf"')));
AVGu = [AASW0,AASW10,AASW30,AASW50];

AASW0 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG mf"')));
AASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv AASW10 mf"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"AVGv AASW30 mf"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGv AASW50 mf"')));
AVGv = [AASW0,AASW10,AASW30,AASW50];

DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
AVG_SPEED = ["(m/s)";"SPEED";SPEED'];

FINAL = [DIST DIRECTION AVG_SPEED SIGM_SPEED];
writematrix(FINAL,[outdir, NewName])
% ---------
% AA TU
NewName =strcat('AASW-AANE_t_TU03A_',ID,'.csv');
% NewName ='AASW-AANE_t_TU03A.csv';
DIST = ["(m)", "DIST", -50, -30, -10,]';

AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW10 t"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW30 t"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED AASW50 t"')));
SPEED = [AASW10,AASW30,AASW50];

AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED AASW10 t"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED AASW30 t"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED AASW50 t"')));
SIGM_SPEED = ["(m/s)","SIGM_SPEED",AASW10,AASW30,AASW50]';

AASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu AASW10 t"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"AVGu AASW30 t"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGu AASW50 t"')));
AVGu = [AASW10,AASW30,AASW50];

AASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv AASW10 t"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"AVGv AASW30 t"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGv AASW50 t"')));
AVGv = [AASW10,AASW30,AASW50];

AASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGw AASW10 t"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"AVGw AASW30 t"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGw AASW50 t"')));
AVGw = [AASW10,AASW30,AASW50];

AASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw AASW10 t"')));
AASW30 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw AASW30 t"')));
AASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw AASW50 t"')));
SIGMw = [AASW10,AASW30,AASW50];
SIGMw = ["(m/s)";"SIGMw";SIGMw'];

DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
AVG_SPEED = ["(m/s)";"SPEED";SPEED'];
HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
AVGw = ["(m/s)";"AVGw"; AVGw'];

FINAL = [DIST DIRECTION AVG_SPEED UPWASH SIGM_SPEED];
writematrix(FINAL,[outdir, NewName])

% Vertical Profile
% ----------------------------------------------------------------
%50m Tower
NewName =strcat('HT_50m_TU03A_vert_',ID,'.csv');
% NewName ='HT_50m_TU03A_vert.csv';
DZ = ["(m)", "DZ", 1, 3, 5, 8, 15, 24, 34, 49]';
HT_1 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 1 m"')));
HT_3 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 3 m"')));
HT_5 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 5 m"')));
HT_8 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 8 m"')));
HT_15 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 15 m"')));
HT_24 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 24 m"')));
HT_34 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 34 m"')));
HT_49 = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT tower 49 m"')));
SPEED = [HT_1,HT_3,HT_5,HT_8,HT_15,HT_24,HT_34,HT_49,];

HT_1 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 1 m"')));
HT_3 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 3 m"')));
HT_5 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 5 m"')));
HT_8 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 8 m"')));
HT_15 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 15 m"')));
HT_24 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 24 m"')));
HT_34 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 34 m"')));
HT_49 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT tower 49 m"')));
SIGM_SPEED = ["(m/s)","SIGM_SPEED",HT_1,HT_3,HT_5,HT_8,HT_15,HT_24,HT_34,HT_49,];

AVG_SPEED = ["(m/s)";"SPEED";SPEED'];

FINAL = [DZ AVG_SPEED SIGM_SPEED'];
writematrix(FINAL,[outdir, NewName])


% ----------------------------------------------------------------
% %CP FRG 17m Tower
% NewName =strcat('CP_FRG17m_TU03A_vert_',ID,'.csv');
% DZ = ["(m)", "DZ", 1.6, 3, 5, 7.4, 10, 16]';
% 
% CP_1 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG 17 m tower 1.6m"')));
% CP_3 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG 17 m tower 3m"')));
% CP_5 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG 17 m tower 5m"')));
% CP_8 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG 17 m tower 7.4m"')));
% CP_10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG 17 m tower 10m"')));
% CP_16 = A.data(slot,find(strcmp(A.colheaders,'"SPEED CP FRG 17 m tower 16m"')));
% SPEED = [CP_1,CP_3,CP_5,CP_8,CP_10,CP_16];
% 
% CP_1 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG 17 m tower 1.6m"')));
% CP_3 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG 17 m tower 3m"')));
% CP_5 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG 17 m tower 5m"')));
% CP_7 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG 17 m tower 7.4m"')));
% CP_10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG 17 m tower 10m"')));
% CP_16 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED CP FRG 17 m tower 16m"')));
% SIGM_SPEED = ["(m/s)", "SIGM SPEED",CP_1,CP_3,CP_5,CP_8,CP_10,CP_16];
% 
% CP_1 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw CP FRG 17 m tower 1.6m"')));
% CP_3 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw CP FRG 17 m tower 3m"')));
% CP_5 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw CP FRG 17 m tower 5m"')));
% CP_8 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw CP FRG 17 m tower 7.4m"')));
% CP_10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw CP FRG 17 m tower 10m"')));
% CP_16 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw CP FRG 17 m tower 16m"')));
% SIGMw = ["(m/s)", "SIGMw", CP_1,CP_3,CP_5,CP_8,CP_10,CP_16];
% 
% CP_1 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG 17 m tower 1.6m"')));
% CP_3 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG 17 m tower 3m"')));
% CP_5 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG 17 m tower 5m"')));
% CP_8 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG 17 m tower 7.4m"')));
% CP_10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG 17 m tower 10m"')));
% CP_16 = A.data(slot,find(strcmp(A.colheaders,'"AVGu CP FRG 17 m tower 16m"')));
% AVGu = [CP_1,CP_3,CP_5,CP_8,CP_10,CP_16];
% 
% CP_1 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG 17 m tower 1.6m"')));
% CP_3 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG 17 m tower 3m"')));
% CP_5 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG 17 m tower 5m"')));
% CP_8 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG 17 m tower 7.4m"')));
% CP_10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG 17 m tower 10m"')));
% CP_16 = A.data(slot,find(strcmp(A.colheaders,'"AVGv CP FRG 17 m tower 16m"')));
% AVGv = [CP_1,CP_3,CP_5,CP_8,CP_10,CP_16];
% 
% CP_1 = A.data(slot,find(strcmp(A.colheaders,'"AVGw CP FRG 17 m tower 1.6m"')));
% CP_3 = A.data(slot,find(strcmp(A.colheaders,'"AVGw CP FRG 17 m tower 3m"')));
% CP_5 = A.data(slot,find(strcmp(A.colheaders,'"AVGw CP FRG 17 m tower 5m"')));
% CP_8 = A.data(slot,find(strcmp(A.colheaders,'"AVGw CP FRG 17 m tower 7.4m"')));
% CP_10 = A.data(slot,find(strcmp(A.colheaders,'"AVGw CP FRG 17 m tower 10m"')));
% CP_16 = A.data(slot,find(strcmp(A.colheaders,'"AVGw CP FRG 17 m tower 16m"')));
% AVGw = [CP_1,CP_3,CP_5,CP_8,CP_10,CP_16];
% 
% % calculate secondary variables and add headers to primary variable columns
% HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
% UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
% DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
% DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
% AVG_SPEED = ["(m/s)";"SPEED";SPEED'];
% 
% % Create final table
% FINAL = [DZ DIRECTION UPWASH AVG_SPEED SIGM_SPEED SIGMw];
% FINAL_TU = FINAL([3 4 6 7],:);
% FINAL_MF = FINAL([5 7 8],[1 2 4]);
% 
% % Write final table to CSV file in OUT directory
% NewName =strcat('CP_FRG17m_mf_TU03A_vert_',ID,'.csv');
% writematrix(FINAL_MF,[outdir, NewName])
% NewName =strcat('CP_FRG17m_t_TU03A_vert_',ID,'.csv');
% writematrix(FINAL_TU,[outdir, NewName])
% 
% % ----------------------------------------------------------------
% %ASW60 vertical
% NewName =strcat('ASW60_UK30m_TU03A_vert_',ID,'.csv');
% DZ = ["(m)", "DZ", 6, 10, 20, 31]';
% 
% UK_6 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW UK 30 m tower 6m"')));
% UK_10= A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW UK 30 m tower 10m"')));
% UK_20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW UK 30 m tower 20m"')));
% UK_31 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW UK 30 m tower 31m"')));
% SPEED = [UK_6, UK_10, UK_20, UK_31];
% 
% UK_6 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW UK 30 m tower 6m"')));
% UK_10= A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW UK 30 m tower 10m"')));
% UK_20 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW UK 30 m tower 20m"')));
% UK_31 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW UK 30 m tower 31m"')));
% SIGM_SPEED = ["(m/s)", "SIGM SPEED", UK_6, UK_10, UK_20, UK_31];
% 
% UK_6 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW UK 30 m tower 6m"')));
% UK_10= A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW UK 30 m tower 10m"')));
% UK_20 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW UK 30 m tower 20m"')));
% UK_31 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW UK 30 m tower 31m"')));
% AVGu = [UK_6, UK_10, UK_20, UK_31];
% 
% UK_6 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW UK 30 m tower 6m"')));
% UK_10= A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW UK 30 m tower 10m"')));
% UK_20 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW UK 30 m tower 20m"')));
% UK_31 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW UK 30 m tower 31m"')));
% AVGv = [UK_6, UK_10, UK_20, UK_31];
% 
% UK_6 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW UK 30 m tower 6m"')));
% UK_10= A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW UK 30 m tower 10m"')));
% UK_20 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW UK 30 m tower 20m"')));
% UK_31 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW UK 30 m tower 31m"')));
% AVGw = [UK_6, UK_10, UK_20, UK_31];
% 
% HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
% UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
% DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
% DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
% AVG_SPEED = ["(m/s)";"SPEED";SPEED'];
% 
% FINAL = [DZ DIRECTION UPWASH AVG_SPEED SIGM_SPEED];
% writematrix(FINAL_MF,[outdir, NewName])
