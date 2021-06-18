% N. Crump  6/11/2021
%Matlabs script to transpose the Askervein hill case files into a form
%usable by the FDS falidation script.

%Mission Process:
%-Read in the decv.csv file from the out directory
%-Transpose the matrix such that each device has it's own row, not column
%-Identify the ASW and ANE (Line A SouthWest -> NorthEast)Devices
%-Insert the Distance from Hiltop column appropriately
%-Sort Rows in Acending order (-100->40)
%-output new devc.csv file

%Questions Checklist:
% -Why did blaze add "..." to the device titles? Is it my fault? it makes
% reading in the column headers very strange -> '"..."' but the code will
% no longer work without them. To check this, >> A.colheaders
% 
% 


outdir = '../../../out/Askervein_Hill/';

%import data from file
filename = 'Askervein_TU03A_16m_devc.csv';
A = importdata([outdir,filename], ",", 2);

%Defines a new name for the output file
NewName = strrep(filename,'devc','devc_read');

%Time_Stop allows the user to change what time they want from the table
%This way only this variable needs to be changed when the cases run longer
TimeStop = 240;
TimeStop1 = TimeStop -1;
TimeStop2 = TimeStop +1;
slot = find( A.data(:,(strcmp(A.colheaders,'Time'))) > TimeStop1 & ...
             A.data(:,(strcmp(A.colheaders,'Time'))) < TimeStop2 );

         
% Create Distance from hiltop column. 
DIST = ["(m)", "DIST", -85, -60, -50, -35, -20, -10, 0, 10, 20, 40]';
         
%data at TIME_STOP (slot variable), this array holds the avarage SPEED data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"SPEED HT Hill Top"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW UK 30 m tower"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SPEED ANE40"')));
SPEED = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

%data at TIME_STOP (slot variable), this array holds the SIGM SPEED data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED HT Hill Top"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW UK 30 m tower"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SIGM SPEED ANE40"')));
% SIGM_SPEED is not used to calculate other values, so header is added here
SIGM_SPEED = ["(m/s)","SIGM SPEED",ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,]';

% %data at TIME_STOP (slot variable), this array holds the SIGMu data
%              -(not currently output by FDS file)-
% A_HT = A.data(slot,find(strcmp(A.colheaders,'"SIGMu HT Hill Top"')));
% ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ASW85"')));
% ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ASW UK 30 m tower"')));
% ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ASW50"')));
% ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ASW35"')));
% ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ASW20"')));
% ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ASW10"')));
% ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ANE10"')));
% ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ANE20"')));
% ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SIGMu ANE40"')));
% SIGMu = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];
% INTu = SIGMu./SPEED;
% SIGMu = ["(m/s)";"SIGMu";SIGMu'];
% INTu = ["NONE";"INTv";INTu'];

% %data at TIME_STOP (slot variable), this array holds the SIGMv data
%              -(not currently output by FDS file)-
% A_HT = A.data(slot,find(strcmp(A.colheaders,'"SIGMv HT Hill Top"')));
% ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ASW85"')));
% ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ASW UK 30 m tower"')));
% ASW50 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ASW50"')));
% ASW35 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ASW35"')));
% ASW20 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ASW20"')));
% ASW10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ASW10"')));
% ANE10 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ANE10"')));
% ANE20 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ANE20"')));
% ANE40 = A.data(slot,find(strcmp(A.colheaders,'"SIGMv ANE40"')));
% SIGMv = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];
% INTv = SIGMv./SPEED;
% SIGMv = ["(m/s)";"SIGMv";SIGMv'];
% INTv = ["NONE";"INTv";INTv'];

%data at TIME_STOP (slot variable), this array holds the SIGMw data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"SIGMw HT Hill Top"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"SIGMw ASW UK 30 m tower"')));
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
INTw = ["NONE";"INTw";INTw'];

%data at TIME_STOP (slot variable), this array holds the AVGu data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"AVGu HT Hill Top"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW UK 30 m tower"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"AVGu ANE40"')));
AVGu = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

%data at TIME_STOP (slot variable), this array holds the AVGv data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"AVGv HT Hill Top"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW UK 30 m tower"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"AVGv ANE40"')));
AVGv = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];

%data at TIME_STOP (slot variable), this array holds the AVGw data
A_HT = A.data(slot,find(strcmp(A.colheaders,'"AVGw HT Hill Top"')));
ASW85 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW85"')));
ASW60 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW UK 30 m tower"')));
ASW50 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW50"')));
ASW35 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW35"')));
ASW20 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW20"')));
ASW10 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ASW10"')));
ANE10 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ANE10"')));
ANE20 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ANE20"')));
ANE40 = A.data(slot,find(strcmp(A.colheaders,'"AVGw ANE40"')));
AVGw = [ASW85,ASW60,ASW50,ASW35,ASW20,ASW10,A_HT,ANE10,ANE20,ANE40,];


% calculate secondary variables and add headers to primary variable columns
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./SPEED)'];
DIRECTION_S0 = atand(AVGv./AVGu).*-1 + 270;
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
AVG_SPEED = ["(m/s)";"SPEED";SPEED'];
AVG_U = ["(m/s)";"AVGu";AVGu'];

% In many graphs there is a column listed as speed. when this column is
% missing there is a column called Ubar. SPEED in FDS is mean speed, but
% I've decided to plot Ubar as well for now

% Create final table
FINAL = [DIST DIRECTION UPWASH AVG_SPEED AVG_U SIGM_SPEED SIGMw INTw];

% Write final table
writematrix(FINAL,[outdir, NewName])


% ----------------------------------------------------------------
% Next up, change file in EXP to add comparison column for SIGM SPEED, and
% to automatically calculate INT values for ASW60

% SIGMspeed = sqrt( SIGMu^2 + SIGMv^2)