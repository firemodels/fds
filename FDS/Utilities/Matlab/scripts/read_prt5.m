% Trettel and McDermott
% 7-5-11
% read_prt5.m
%
% This function reads the FDS 'part' file (*.prt5).
%
% Example:
%
% >> read_prt5('terminal_velocity.prt5','real*8')
%
% precision = 'real*4' if EB_PART_FILE=.FALSE. on the DUMP line
% precision = 'real*8' if EB_PART_FILE=.TRUE. on the DUMP line

function [STIME, XP, YP, ZP, QP] = read_prt5(filename,precision)

fid = fopen(filename);

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

for NPC=1:N_PART
    
    DUMMY = fread(fid,1,'integer*4');
    PC    = fread(fid,2,'integer*4'); N_QUANTITIES(NPC) = PC(1);
    DUMMY = fread(fid,1,'integer*4');
    
    for NQ=1:N_QUANTITIES(NPC)
        DUMMY               = fread(fid,1,'integer*4');
        SMOKEVIEW_LABEL{NQ} = fgets(fid,30);
        DUMMY               = fread(fid,1,'integer*4');
        
        DUMMY     = fread(fid,1,'integer*4');
        UNITS{NQ} = fgets(fid,30);
        DUMMY     = fread(fid,1,'integer*4');
    end
    
end

n = 0;
while ~feof(fid)
    n = n + 1;
    
    DUMMY = fread(fid,1,'integer*4');
    stime_tmp = fread(fid,1,precision);
    DUMMY = fread(fid,1,'integer*4');
    
    if size(stime_tmp,1)==0
        break
    else
        STIME(n) = stime_tmp;
    end
    
    for NPC=1:N_PART
        
        DUMMY = fread(fid,1,'integer*4');
        NPLIM = fread(fid,1,'integer*4');
        DUMMY = fread(fid,1,'integer*4');
        
        DUMMY = fread(fid,1,'integer*4');
        xp = fread(fid,NPLIM,precision);
        yp = fread(fid,NPLIM,precision);
        zp = fread(fid,NPLIM,precision);
        DUMMY = fread(fid,1,'integer*4');
        
        for NP=1:NPLIM
            XP(n,NP,NPC) = xp(NP);
            YP(n,NP,NPC) = yp(NP);
            ZP(n,NP,NPC) = zp(NP);
        end
        %clear xp yp zp
        
        DUMMY = fread(fid,1,'integer*4');
        TA    = fread(fid,NPLIM,'integer*4');
        DUMMY = fread(fid,1,'integer*4');
        
        if N_QUANTITIES(NPC)>0
            DUMMY = fread(fid,1,'integer*4');
            for NQ=1:N_QUANTITIES(NPC)
                qp = fread(fid,NPLIM,precision);
                for NP=1:NPLIM
                    QP(n,NP,NPC,NQ) = qp(NP);
                end
                %clear qp
            end
            DUMMY = fread(fid,1,'integer*4');
        end
        
    end
    
end
fclose(fid);

%display('Part file read successfully!')

% Examples for plotting position and quantities

%plot(STIME,ZP(:,1,1))   % XP(time step range, particle #, part class #)
%size(STIME)
%size(QP(:, 1, 1, 1))
%plot(STIME,QP(:,1,1,1)) % QP(time step range, particle #, part class #, quantity #)
%min(QP(:,1,1,1))
