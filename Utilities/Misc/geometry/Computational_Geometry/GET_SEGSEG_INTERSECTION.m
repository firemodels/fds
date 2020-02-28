function [SVARV,SLENV,INT_FLG]=GET_SEGSEG_INTERSECTION(P1,D1,P2,D2)

global IAXIS JAXIS NOD1 NOD2 EDG1 EDG2 GEOMEPS

SVARV(NOD1:NOD2,EDG1:EDG2) = -1.;
SLENV(EDG1:EDG2)           = -1.;

% Test for segment-segment intersection:
E(IAXIS:JAXIS) = P2(IAXIS:JAXIS) - P1(IAXIS:JAXIS);
KRS = D1(IAXIS)*D2(JAXIS) - D1(JAXIS)*D2(IAXIS); KRS2=KRS^2.;
L12 = D1(IAXIS)^2. + D1(JAXIS)^2.;
L22 = D2(IAXIS)^2. + D2(JAXIS)^2.;
         
% Case of segments not parallel.
if ( KRS2 > GEOMEPS^2.*L12*L22)
    SVR = (E(IAXIS)*D2(JAXIS)-E(JAXIS)*D2(IAXIS))/ KRS;
    if(SVR < -GEOMEPS || (SVR-1.) > GEOMEPS)
        % intersection not a point of segment SEG.
        INT_FLG=0;
        return
    end
    TVR = (E(IAXIS)*D1(JAXIS)-E(JAXIS)*D1(IAXIS))/ KRS;
    if(TVR < -GEOMEPS || (TVR-1.) > GEOMEPS)
        % intersection not a point of segment SEG2.
        INT_FLG=0;
        return
    end

    % Intersection a point on SEG and SEG2.
    SLENV(EDG1)      = sqrt(L12);
    SLENV(EDG2)      = sqrt(L22);
    SVARV(NOD1,EDG1) = SVR*SLENV(EDG1);
    SVARV(NOD1,EDG2) = TVR*SLENV(EDG2);
    INT_FLG=1;
    return
    
end

% Parallel Segments:
E2 = E(IAXIS)^2. + E(JAXIS)^2.;
KRS= E(IAXIS)*D1(JAXIS) - E(JAXIS)*D1(IAXIS); KRS2=KRS^2.;
if( KRS2 > GEOMEPS^2.*L12*E2 )
    % Segments are different.
    INT_FLG=0;
    return
end
% Segment lines are the same. Overlap tests.
S1  = dot(D1,E)/L12; S2 = S1+dot(D1,D2)/L12;
SMIN=min(S1,S2); SMAX=max(S1,S2);
if( (1.+GEOMEPS) < SMIN || (0.-GEOMEPS) > SMAX)
    INT_FLG=0;
    return
end

SLENV(EDG1)      = sqrt(L12);
SLENV(EDG2)      = sqrt(L22);

if( (1.+GEOMEPS) > SMIN ) % SMIN between P1 and P1+D1
    if ( (0.-GEOMEPS) < SMAX) % SMAX greater that P1
       if (0. < SMIN) % SMIN higher that P1
           SVARV(NOD1,EDG1) = SMIN*SLENV(EDG1); % First crossing on P1-P1+D1
           if (abs(SMIN-S1) < GEOMEPS/2.) % SMIN is P2
               SVARV(NOD1,EDG2)=0.; % First crossing in P2-P2+D2
           else % SMIN is P2+D2
               SVARV(NOD2,EDG2)=1.*SLENV(EDG2); % Second crossing in P2-P2+D2
           end
       else % SMIN lower than P1
           SVARV(NOD1,EDG1) = 0.; % First crossing in P1-P1+D1
           if (abs(SMIN-S1) < GEOMEPS/2.) % SMIN os P2
               SVARV(NOD1,EDG2)=-SMIN*SLENV(EDG1); % First crossing in P2-P2-D2
           else
               SVARV(NOD2,EDG2)=SMAX*SLENV(EDG1);
           end
       end
       if (1. > SMAX)
           SVARV(NOD2,EDG1) = SMAX*SLENV(EDG1);
           if(abs(SMAX-S1) < GEOMEPS/2.) % SMAX is P2
               SVARV(NOD1,EDG2)=0.*SLENV(EDG2);
           else
               SVARV(NOD2,EDG2)=1.*SLENV(EDG2);
           end
       else
           SVARV(NOD2,EDG1) = 1.*SLENV(EDG1);
           if(abs(SMAX-S1) < GEOMEPS/2.) % SMAX is P2
               SVARV(NOD1,EDG2)=(SMAX-1.)*SLENV(EDG1);
           else
               SVARV(NOD2,EDG2)=(1.-SMIN)*SLENV(EDG1);
               
           end
       end
       INT_FLG=2;
       return
    else
       % SMAX = 0.
       SVARV(NOD1,EDG1) = 0.;
       if(abs(SMAX-S1) < GEOMEPS/2.)
           SVARV(NOD1,EDG2) = 0.;
       else
           SVARV(NOD1,EDG2) = 1.*SLENV(EDG2);
       end 
       INT_FLG=1;
       return     
    end
else
    % SMIN = 1.
    SVARV(NOD1,EDG1) = 1.*SLENV(EDG1);
    if(abs(SMIN-S1) < GEOMEPS/2.)
        SVARV(NOD1,EDG2) = 0.;
    else
        SVARV(NOD1,EDG2) = 1.*SLENV(EDG2);
    end 
    INT_FLG=1;
    return
end

return