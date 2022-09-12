% N. Crump  6/11/2021 - 8/5/2021
%Matlabs script to transpose the Askervein hill case files into a form
%usable by the FDS validation script.

% Script Process:
%-Read in the decv.csv file from the out directory
%-Build the output matrix with sorted rows
%-Insert the Distance from Hiltop column
%-output new .csv file

% Navigation
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

% Questions Checklist:
% reading in the column headers is very strange '"..."' the code will
% not work without the extra "". To see this, >> A.colheaders
 
% Debugging notes:
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

%Time_Stop, set up to  take the latest value on the table
Times = A.data(:,(strcmp(A.colheaders,'Time')));
TimeStop = max(Times);
TimeStop1 = TimeStop - 1;
TimeStop2 = TimeStop + 1;
slot = find( A.data(:,(strcmp(A.colheaders,'Time'))) > TimeStop1 & ...
             A.data(:,(strcmp(A.colheaders,'Time'))) < TimeStop2 );

% A Line 10m height
% ----------------------------------------------------------------
NewName = strcat('ASW-ANE_TU03A_',ID,'.csv');
% Define Towers and Quantities
DIST = ["(m)", "DIST", -85, -60, -50, -35, -20, -10, 0, 10, 20, 40]';
TOWERS = ["ASW85","ASW UK 30 m tower 10m","ASW50","ASW35","ASW20","ASW10","HT 10 m t","ANE10","ANE20","ANE40"];
QUANTITY = ["SPEED ","SIGM SPEED ","SIGMw ","AVGu ","AVGv ","AVGw "];

colls = [];
f = '"';
for this_quant = QUANTITY
    new_row = [];
    for this_tower = TOWERS
        this_name = (append(f,this_quant,this_tower,f));
        value = A.data(slot,find(strcmp(A.colheaders,this_name)));
        new_row = [new_row,value];
    end
    colls = cat(1,colls,new_row);
end

AVG_SPEED = ["(m/s)";"SPEED";colls(1,:)'];
SIGM_SPEED = ["(m/s)";"SIGM_SPEED"; colls(2,:)'];
SIGMw = ["(m/s)";"SIGMw";colls(3,:)'];
AVGu = colls(4,:);
AVGv = colls(5,:);
AVGw = colls(6,:);

% calculate secondary variables and add headers to primary variable columns
HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];


% Create final table
FINAL = [DIST DIRECTION UPWASH AVG_SPEED SIGM_SPEED SIGMw];
% Write final table to CSV file in OUT directory
writematrix(FINAL,[outdir, NewName])

% AA Line
% ----------------------------------------------------------------
%SPEED
NewName = strcat('AASW-AANE_TU03A_',ID,'.csv');
%NewName = 'AASW-AANE_TU03A_10m.csv';
DIST = ["(m)","DIST",-90,-80,-70,-60,-50,-40,-30,-20,-10,0,0,10,20,30,40,50,60]';
TOWERS = ["AASW90","AASW80","AASW70","AASW60","AASW50 mf","AASW40","AASW30 mf","AASW20","AASW10 mf"];
TOWERS = [TOWERS,"CP BSE 40","CP UK mf","AANE10","AANE20","AANE30","AANE40","AANE50","AANE60"];

    SPEED = [];
for this_tower = TOWERS
    this_name = convertStringsToChars(append(f,"SPEED ",this_tower,f));
    value = A.data(slot,find(strcmp(A.colheaders,this_name)));
    SPEED = [SPEED,value];
end

FINAL = [DIST ["(m/s)";"SPEED";SPEED']];
writematrix(FINAL,[outdir, NewName])

% B Line
% ----------------------------------------------------------------
%SPEED
NewName = strcat('BNW-BSE_TU03A_',ID,'.csv');
% NewName = 'BNW-BSE_TU03A_10m.csv';
DIST = ["(m)","DIST",-20,-10,0,10,20,30,40,40,50,60,70,80,90,100,110,130]';
TOWERS = ["BNW20","BNW10","HT 10 m mf","BSE10","BSE20","BSE30","BSE40","CP UK mf"];
TOWERS = [TOWERS,"BSE50","BSE60","BSE70","BSE80","BSE90","BSE100","BSE110","BSE130"];

    SPEED = [];
for this_tower = TOWERS
    this_name = append(f,"SPEED ",this_tower,f);
    value = A.data(slot,find(strcmp(A.colheaders,this_name)));
    SPEED = [SPEED,value];
end

AVG_SPEED = ["(m/s)";"SPEED";SPEED'];
FINAL = [DIST AVG_SPEED];
writematrix(FINAL,[outdir, NewName])

% AA Short Line mf and tu Mean Flow
% ----------------------------------------------------------------
%AA MF
NewName = strcat('AASW-AANE_mf_TU03A_',ID,'.csv');
% NewName = 'AASW-AANE_mf_TU03A.csv';
DIST = ["(m)", "DIST", -50, -30, -10, 0]';
TOWERS = ["CP FRG mf","AASW10 mf","AASW30 mf","AASW50 mf"];
QUANTITY = ["SPEED ","SIGM SPEED ","AVGu ","AVGv "];

colls = [];
for this_quant = QUANTITY
    new_row = [];
    for this_tower = TOWERS
        this_name = append(f,this_quant,this_tower,f);
        value = A.data(slot,find(strcmp(A.colheaders,this_name)));
        new_row = [new_row,value];
    end
    colls = cat(1,colls,new_row);
end

AVG_SPEED = ["(m/s)";"SPEED";colls(1,:)'];
SIGM_SPEED = ["(m/s)","SIGM_SPEED",colls(2,:)]';
AVGu = colls(3,:);
AVGv = colls(4,:);

DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];

FINAL = [DIST DIRECTION AVG_SPEED SIGM_SPEED];
writematrix(FINAL,[outdir, NewName])

% ---------
% AA TU
NewName = strcat('AASW-AANE_t_TU03A_',ID,'.csv');
% NewName = 'AASW-AANE_t_TU03A.csv';
DIST = ["(m)", "DIST", -50, -30, -10,]';
TOWERS = ["AASW10 t","AASW30 t","AASW50 t"];
QUANTITY = ["SPEED ","SIGM SPEED ","SIGMw ","AVGu ","AVGv ","AVGw "];

colls = [];
for this_quant = QUANTITY
    new_row = [];
    for this_tower = TOWERS
        this_name = append(f,this_quant,this_tower,f);
        value = A.data(slot,find(strcmp(A.colheaders,this_name)));
        new_row = [new_row,value];
    end
    colls = cat(1,colls,new_row);
end

AVG_SPEED = ["(m/s)";"SPEED";colls(1,:)'];
SIGM_SPEED = ["(m/s)";"SIGM SPEED";colls(2,:)'];
SIGMw = ["(m/s)";"SIGMw";colls(3,:)'];
AVGu = colls(4,:);
AVGv = colls(5,:);
AVGw = colls(6,:);

DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];
HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
AVGw = ["(m/s)";"AVGw"; AVGw'];

FINAL = [DIST DIRECTION AVG_SPEED UPWASH SIGM_SPEED];
writematrix(FINAL,[outdir, NewName])

% Vertical Profile
% ----------------------------------------------------------------
%50m Tower
NewName = strcat('HT_50m_TU03A_vert_',ID,'.csv');
% NewName = 'HT_50m_TU03A_vert.csv';
DZ = ["(m)", "DZ", 1, 3, 5, 8, 15, 24, 34, 49]';
TOWERS = ["1 m","3 m","5 m","8 m","15 m","24 m","34 m","49 m"];
QUANTITY = ["SPEED ","SIGM SPEED "];

colls = [];
for this_quant = QUANTITY
    new_row = [];
    for this_tower = TOWERS
        this_name = append(f,this_quant,"HT tower ",this_tower,f);
        value = A.data(slot,find(strcmp(A.colheaders,this_name)));
        new_row = [new_row,value];
    end
    colls = cat(1,colls,new_row);
end

AVG_SPEED = ["(m/s)";"SPEED";colls(1,:)'];
SIGM_SPEED = ["(m/s)";"SIGM_SPEED";colls(2,:)'];

FINAL = [DZ AVG_SPEED SIGM_SPEED];
writematrix(FINAL,[outdir, NewName])
%comment below
% ----------------------------------------------------------------
%CP FRG 17m Tower
DZ = ["(m)", "DZ", 1.6, 3, 5, 7.4, 10, 16]';
TOWERS = ["1.6m","3m","5m","7.4m","10m","16m"];
QUANTITY = ["SPEED ","SIGM SPEED ","SIGMw ","AVGu ","AVGv ","AVGw "];

colls = [];
for this_quant = QUANTITY
    new_row = [];
    for this_tower = TOWERS
        this_name = append(f,this_quant,"CP FRG 17 m tower ",this_tower,f);
        value = A.data(slot,find(strcmp(A.colheaders,this_name)));
        new_row = [new_row,value];
    end
    colls = cat(1,colls,new_row);
end

AVG_SPEED = ["(m/s)";"SPEED";colls(1,:)'];
SIGM_SPEED = ["(m/s)", "SIGM SPEED",colls(2,:)]';
SIGMw = ["(m/s)", "SIGMw",colls(3,:)]';
AVGu = colls(4,:);
AVGv = colls(5,:);
AVGw = colls(6,:);

% calculate secondary variables and add headers to primary variable columns
HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];

% Create final tables
FINAL = [DZ DIRECTION UPWASH AVG_SPEED SIGM_SPEED SIGMw];
FINAL_TU = FINAL([3 4 6 7],:);
FINAL_MF = FINAL([5 7 8],[1 2 4]);

% Write final table to CSV file in OUT directory
NewName = strcat('CP_FRG17m_mf_TU03A_vert_',ID,'.csv');
writematrix(FINAL_MF,[outdir, NewName])
NewName = strcat('CP_FRG17m_t_TU03A_vert_',ID,'.csv');
writematrix(FINAL_TU,[outdir, NewName])

% ----------------------------------------------------------------
%ASW60 vertical
NewName = strcat('ASW60_UK30m_TU03A_vert_',ID,'.csv');
DZ = ["(m)", "DZ", 6, 10, 20, 31]';
TOWERS = ["6m","10m","20m","31m"];
QUANTITY = ["SPEED ","SIGM SPEED ","AVGu ","AVGv ","AVGw "];

colls = [];
for this_quant = QUANTITY
    new_row = [];
    for this_tower = TOWERS
        this_name = append(f,this_quant,"ASW UK 30 m tower ",this_tower,f);
        value = A.data(slot,find(strcmp(A.colheaders,this_name)));
        new_row = [new_row,value];
    end
    colls = cat(1,colls,new_row);
end

AVG_SPEED = ["(m/s)";"SPEED";colls(1,:)'];
SIGM_SPEED = ["(m/s)", "SIGM SPEED", colls(2,:)]';
AVGu = colls(3,:);
AVGv = colls(4,:);
AVGw = colls(5,:);

HORIZ_VELO = (AVGv.^2 + AVGu.^2).^.5;
UPWASH = ["(deg)";"UPWASH"; atand(AVGw./HORIZ_VELO)'];
DIRECTION_S0 = (AVGu./abs(AVGu)).*90 + 180 - atand(AVGv./AVGu); 
DIRECTION = ["(deg)";"DIRECTION"; DIRECTION_S0'];

FINAL = [DZ DIRECTION UPWASH AVG_SPEED SIGM_SPEED];
writematrix(FINAL,[outdir, NewName])