function [IFLG,NPOLY,ILO_POLY,NSG_POLY,NSEG,SEG_CELL,SEG_POS]=GET_CLOSED_POLYLINES(NSEG,SEG_CELL,SEG_POS)

global NOD1 NOD2 IBM_MAX_WSTRIANG_SEG

IFLG      =  1;
MAX_POLY  = 20;
SEG_CELL2 = zeros(NOD2+IBM_MAX_WSTRIANG_SEG+2,NSEG);
SEG_POS2  = zeros(1,NSEG);
COUNTED   = zeros(1,NSEG);

% First collapse segments to most frequent body:
NBOD=1;
BOD(NBOD)=SEG_CELL(6,1);
for ISEG=2:NSEG
    INLIST=false;
    for IBOD=1:NBOD
        if(SEG_CELL(6,ISEG) == BOD(IBOD))
            INLIST=true;
            break
        end
    end
    if(~INLIST)
        NBOD=NBOD+1;
        BOD(NBOD)=SEG_CELL(6,ISEG);
    end
end
CTBOD = zeros(1,NBOD);
for IBOD=1:NBOD
    for ISEG=1:NSEG
        if(SEG_CELL(6,ISEG) == BOD(IBOD))
            CTBOD(IBOD) = CTBOD(IBOD) + 1;
        end
    end
end
[MXBOD,MIBOD]=max(CTBOD);

for ISEG=1:NSEG
   if(COUNTED(ISEG)~=0); continue; end
   CISEG=0;
   for ISEG2=1:NSEG
       if (COUNTED(ISEG2)~=0); continue; end
       if( ISEG2==ISEG ); continue; end
       if( SEG_CELL(NOD1,ISEG) == SEG_CELL(NOD1,ISEG2) && ...
           SEG_CELL(NOD2,ISEG) == SEG_CELL(NOD2,ISEG2) )
          if(SEG_CELL(6,ISEG) == BOD(MIBOD))
              % ISEG should be COUNTED +1; ISEG2 -1.
              COUNTED(ISEG) = 1;
              COUNTED(ISEG2)=-1;
              CISEG=1;
          else
              % ISEG should be COUNTED -1; ISEG2 +1.
              COUNTED(ISEG) =-1;
              COUNTED(ISEG2)= 1;
              CISEG=1;
          end
       end
   end
   if(CISEG==0); COUNTED(ISEG) = 1; end
end
NEWSEG=0;
for ISEG=1:NSEG
    if(COUNTED(ISEG)~=1); continue; end
    NEWSEG = NEWSEG +1;
    SEG_CELL2(:,NEWSEG) = SEG_CELL(:,ISEG);
    SEG_POS2(NEWSEG)    = SEG_POS(ISEG);
end
NSEG      = NEWSEG;
SEG_CELL  = SEG_CELL2;
SEG_POS   = SEG_POS2;


% Now make closed polylines:
SEG_CELL2 = zeros(NOD2+IBM_MAX_WSTRIANG_SEG+2,NSEG);
SEG_POS2  = zeros(1,NSEG);
COUNTED   = zeros(1,NSEG);
NPOLY     = 0;
SEG_POLY  = zeros(1,NSEG); % Polyline number for the segment.
ILO_POLY  = zeros(1,MAX_POLY);
NSG_POLY  = zeros(1,MAX_POLY); 
SEG_LEFT  = NSEG;
while(1) % This exterior while loop defined closed polylines in the cell.
    
    % Count one more polyline:
    NPOLY = NPOLY + 1;
    if(NPOLY==1)
        ILO_POLY(NPOLY) = 0;
    else
        ILO_POLY(NPOLY) = ILO_POLY(NPOLY-1) + NSG_POLY(NPOLY-1);
    end
    
    % Find first segment of next polyline:
    FOUNDSEG= 0;
    for ISEG=1:NSEG
        if(~COUNTED(ISEG))
            FOUNDSEG= 1; 
            break
        end
    end
    if(~FOUNDSEG); break; end % Escape if there are no new segments.
    
    % Create new closed polyline:
    NEWSEG   = ILO_POLY(NPOLY)+1;
    SEG_CELL2(:,NEWSEG) = SEG_CELL(:,ISEG);
    SEG_POS2(NEWSEG) = SEG_POS(ISEG);
    COUNTED(ISEG) = true;
    STNOD  = SEG_CELL2(NOD1,NEWSEG);
    PIVNOD = SEG_CELL2(NOD2,NEWSEG); % Pivot Vertex, used to find next segment.
    NSG_POLY(NPOLY) =  NSG_POLY(NPOLY) + 1;
    SEG_POLY(NEWSEG)  = NPOLY;
    SEG_LEFT = SEG_LEFT - 1;
    for NEWSEG = ILO_POLY(NPOLY)+2:NSEG
        FOUNDSEG= 0;
        for ISEG=1:NSEG
            
            if(COUNTED(ISEG)); continue; end

            if(SEG_CELL(NOD1,ISEG)==PIVNOD) % Found the next segment
                FOUNDSEG = 1;
                SEG_CELL2(:,NEWSEG) = SEG_CELL(:,ISEG);
                SEG_POS2(NEWSEG) = SEG_POS(ISEG);
                COUNTED(ISEG) = true;
                PIVNOD = SEG_CELL2(NOD2,NEWSEG); % Pivot Vertex, used to find next segment.
                NSG_POLY(NPOLY) =  NSG_POLY(NPOLY) + 1;
                SEG_POLY(NEWSEG)  = NPOLY;
                SEG_LEFT = SEG_LEFT - 1;
                break
            elseif(SEG_CELL(NOD2,ISEG)==PIVNOD) % Found the next segment
                FOUNDSEG = 1;
                SEG_CELL2(:,NEWSEG) = [SEG_CELL(NOD2,ISEG) SEG_CELL(NOD1,ISEG) SEG_CELL(3:6,ISEG)']';
                SEG_POS2(NEWSEG) = SEG_POS(ISEG);
                COUNTED(ISEG) = true;
                PIVNOD = SEG_CELL2(NOD2,NEWSEG); % Pivot Vertex, used to find next segment.
                NSG_POLY(NPOLY) =  NSG_POLY(NPOLY) + 1;
                SEG_POLY(NEWSEG)  = NPOLY;
                SEG_LEFT = SEG_LEFT - 1;
                break
            end
            
        end
        % Check if for this NEWSEG we didn't find an ISEG:
        if (~FOUNDSEG)
            break
        end
    end
    % Finally, test if polyline is closed:
    if( SEG_CELL2(NOD2,ILO_POLY(NPOLY)+NSG_POLY(NPOLY)) ~= STNOD )
        disp('Polyline not closed')
        return % Return with IFLG = 1.
    end
    % End of new polyline creation.
    % Here if we have less that 3 segments not counted exit while loop.
    if(SEG_LEFT < 3); break; end
    
end

% Per polyline, move last SEG if SEG-1 is different body number:
for IPOLY=1:NPOLY
    FOUND_CHG=0;
    ILO      =ILO_POLY(IPOLY)+1;
    IHI      =ILO_POLY(IPOLY)+NSG_POLY(IPOLY);
    CT       =0;
    for ISEG=ILO:IHI-1
        CT=CT+1;
        if(SEG_CELL2(6,ISEG) ~= SEG_CELL2(6,ISEG+1))
            FOUND_CHG=1;
            break
        end
    end
    if(FOUND_CHG)
        SEG_CELL(:,ILO:IHI-CT)      = SEG_CELL2(:,ISEG+1:IHI);
        SEG_POS(ILO:IHI-CT)         =  SEG_POS2(1,ISEG+1:IHI);
        SEG_CELL(:,IHI-CT+1:IHI) = SEG_CELL2(:,ILO:ISEG);
        SEG_POS(IHI-CT+1:IHI)    =  SEG_POS2(1,ILO:ISEG);        
    else
        SEG_CELL(:,ILO:IHI)= SEG_CELL2(:,ILO:IHI);
        SEG_POS(ILO:IHI)   = SEG_POS2(1,ILO:IHI);
    end
end

% Finally cycle segments to redefine polylines (case of two or more polys
% sharing one point.
STNOD=SEG_CELL(NOD1,1);
NPOLY=1;
COUNT=1;
for ISEG=2:NSEG
   COUNT=COUNT+1;
   SEG_POLY(ISEG)=NPOLY;
   if(SEG_CELL(NOD2,ISEG)==STNOD)
       NSG_POLY(NPOLY) = COUNT;
       if(ISEG==NSEG); break; end
       NPOLY=NPOLY+1;
       ILO_POLY(NPOLY) = ILO_POLY(NPOLY-1) + NSG_POLY(NPOLY-1);
       COUNT=0; STNOD=SEG_CELL(NOD1,ISEG+1);
   end
end

IFLG=0;

return