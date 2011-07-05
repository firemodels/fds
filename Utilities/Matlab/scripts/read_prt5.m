% McDermott and Trettel
% 7-5-11
% read_prt5.m
%
% This script reads the FDS 'part' file (*.prt5).
%
% Currently, only works for reading the particle time and position.
% This script will not work if any QUANTITIES are present in the part file.
% The capability to read the QUANTITIES is coming soon.
%
% Usage: Copy this script to your working directory and modify the filename
% on the fopen line.

close all
clear all

fid = fopen('terminal_velocity.prt5'); % modify this line

% The DUMMY lines are 4 byte placeholders that apparently fortran puts at the
% beginning and end of all lines.  I only knew of this thanks to Glenn.

DUMMY       = fread(fid,1,'integer*4');
ONE_INTEGER = fread(fid,1,'integer*4');
DUMMY       = fread(fid,1,'integer*4');

DUMMY       = fread(fid,1,'integer*4');
INT_VERSION = fread(fid,1,'integer*4');
DUMMY       = fread(fid,1,'integer*4');

DUMMY  = fread(fid,1,'integer*4');
N_PART = fread(fid,1,'integer*4');
DUMMY  = fread(fid,1,'integer*4');

for NP=1:N_PART
    
    DUMMY = fread(fid,1,'integer*4');
    PC    = fread(fid,2,'integer*4'); N_QUANTITIES(NP) = PC(1);
    DUMMY = fread(fid,1,'integer*4');
    
    for NN=1:N_QUANTITIES(NP)
        DUMMY               = fread(fid,1,'integer*4');
        SMOKEVIEW_LABEL{NP} = fgets(fid,30);
        DUMMY               = fread(fid,1,'integer*4');
        
        DUMMY     = fread(fid,1,'integer*4');
        UNITS{NP} = fgets(fid,30);
        DUMMY     = fread(fid,1,'integer*4');
    end
    
end

n = 0;
while ~feof(fid)
    n = n + 1;
    
    DUMMY = fread(fid,1,'integer*4');
    stime_tmp = fread(fid,1,'real*4');
    DUMMY = fread(fid,1,'integer*4');
    
    if size(stime_tmp,1)==0
        break
    else
        STIME(n) = stime_tmp;
    end
    
    for NP=1:N_PART
        
        DUMMY = fread(fid,1,'integer*4');
        NPLIM = fread(fid,1,'integer*4');
        DUMMY = fread(fid,1,'integer*4');
        
        DUMMY = fread(fid,1,'integer*4');
        XP(n,NP) = fread(fid,NPLIM,'real*4');
        YP(n,NP) = fread(fid,NPLIM,'real*4');
        ZP(n,NP) = fread(fid,NPLIM,'real*4');
        DUMMY = fread(fid,1,'integer*4');
        
        DUMMY = fread(fid,1,'integer*4');
        TA    = fread(fid,NPLIM,'integer*4');
        DUMMY = fread(fid,1,'integer*4');
        
        if N_QUANTITIES(NP)>0
            DUMMY = fread(fid,1,'integer*4');
            for NN=1:N_QUANTITIES(NP)
                QP(1:NPLIM,NN,NP) = fread(fid,NPLIM,'real*4');
            end
            DUMMY = fread(fid,1,'integer*4');
        end
        
    end
    
end
fclose(fid);

display('Part file read successfully!')

%plot(STIME,ZP(:,1))
%plot(STIME,QP(:,1,1))
