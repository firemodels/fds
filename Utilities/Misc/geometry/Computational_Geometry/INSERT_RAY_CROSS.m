function [ierr]=INSERT_RAY_CROSS(SVARI,ICRSI,SCRSI,STANI)

global IAXIS JAXIS IBM_SOLID IBM_GASPHASE
global IBM_N_CRS 
global LOW_IND HIGH_IND
global IBM_SVAR_CRS IBM_IS_CRS2 IBM_SEG_CRS IBM_SEG_TAN IBM_BDNUM_CRS

ierr = 1;

IBM_N_CRS = IBM_N_CRS + 1;

% Add in place, ascending value order:
for ICRS=1:IBM_N_CRS % The updated IBM_N_CRS is for ICRS to reach the
                     % initialization value IBM_SVAR_CRS(ICRS)=1/GEOMEPS.
   if ( SVARI < IBM_SVAR_CRS(ICRS) ); break; end
end

% Here copy from the back (updated IBM_N_CRS) to the ICRS location:
% if ICRS=IBM_N_CRS -> nothing gets copied:
for IDUM = IBM_N_CRS : -1 : ICRS+1
   IBM_SVAR_CRS(IDUM)           = IBM_SVAR_CRS(IDUM-1);
   IBM_IS_CRS2(LOW_IND:HIGH_IND+1,IDUM)   = IBM_IS_CRS2(LOW_IND:HIGH_IND+1,IDUM-1);
   IBM_SEG_CRS(IDUM)            = IBM_SEG_CRS(IDUM-1);
   IBM_SEG_TAN(IAXIS:JAXIS,IDUM)= IBM_SEG_TAN(IAXIS:JAXIS,IDUM-1);
   IBM_BDNUM_CRS(IDUM) = IBM_BDNUM_CRS(IDUM-1);
end

IBM_SVAR_CRS(ICRS)             = SVARI;              % x2 location.
IBM_IS_CRS2(LOW_IND:HIGH_IND+1,ICRS)     = ICRSI(LOW_IND:HIGH_IND+1);    % Does point separate GASPHASE from SOLID?
IBM_SEG_CRS(ICRS)              = SCRSI;              % Segment on BOINT_PLANE the crossing belongs to.
IBM_SEG_TAN(IAXIS:JAXIS,ICRS)  = STANI(IAXIS:JAXIS); % IBM_SEG_TAN might not be needed in new implementation.
IBM_BDNUM_CRS(ICRS) = 0;
if (SCRSI > 0)
   if (ICRSI(LOW_IND) == IBM_GASPHASE && ICRSI(HIGH_IND) == IBM_SOLID)
      IBM_BDNUM_CRS(ICRS)  = 1;
   elseif(ICRSI(LOW_IND) == IBM_SOLID && ICRSI(HIGH_IND) == IBM_GASPHASE)
      IBM_BDNUM_CRS(ICRS) =-1;
   end
end

ierr=0;

return


